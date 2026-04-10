##### Script to reproject a batch of homolosine GEOTiffs from ISRIC using gdalwarp #####
##### Camila Medeiros #####
##### 7th Feb 2023 #####
# https://www.isric.org/explore/soilgrids/faq-soilgrids
# https://www.isric.org/explore/soilgrids/faq-soilgrids-2017
# https://stackoverflow.com/questions/55409398/reprojecting-list-of-raster-stacks

## Load packages:
library(raster)
library(rgdal)
library(gdalUtils)


## Set WD:
setwd("/Volumes/CamilaBckup/climate/")

## Load environment layers:
# get path to all GEOTiff files in the specified folder
soil_list <- list.files("ISRICsoilgrids_2023", pattern=".tif$", full.names = F) # set full.names to F to facilitate the use of this list for the loop; would set as T otherwise

# Now we have a list of 61 soil rasters in Homolosine projection at 1 km resolution (.5 degrees).
# However, our species occurrence shapefiles and all other climate variables are projected as 
# Spherical Mercator (EPSG 4326), so we need to re-project them. There are two main ways to go about doing this, one using 
# raster::projectRaster and another using rgdal::gdalwarp. The latter is much faster and yields a product with lower
# error/differences from the original raster (https://gis.stackexchange.com/questions/154835/difference-between-gdalwarp-and-projectraster).
# Thus, here I decided to use gdalwarp to reproject. -- Will need to go back soon to replace functions with the equivalent from sf and terra (start here: https://r-tmap.github.io/tmap-book/geodata.html#crs)

## Define the folder where the reprojected GEOTiffs should be saved in (important: write the *full* file paths here!!):
indir <- "/Volumes/CamilaBckup/climate/ISRICsoilgrids_2023/"
outdir <- "/Volumes/CamilaBckup/climate/ISRICsoilgrids_2023_reprojected/"

## Define the path to the GEOTiff I want to reproject (write full file path here as well)
inputRaster <- paste0(indir, soil_list[60])    # written for loop 
outputRaster <- paste0(outdir, "r", soil_list[60])   # written for loop

## Now I can reproject a single GEOTiff:
gdalwarp(inputRaster, dstfile= outputRaster,
         t_srs= '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', #target crs: EPSG 4326 declaration (Spherical Mercator)
         r ="bilinear", # recommended method to resample continuous rasters, such as ours
         output_Raster=TRUE,
         overwrite=TRUE, 
         verbose=TRUE)

## And here I can reproject all the GEOTiffs in the input directory (might take a while!):
nn <- length(soil_list) # number of GEOTiffs in input directory

for(i in 1:nn){
  inputRaster[i] <- paste0(indir, soil_list[i])   
  
  outputRaster[i] <- paste0(outdir, "r", soil_list[i])
  
  gdalwarp(inputRaster[i], dstfile= outputRaster[i],
           t_srs= '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs', #target crs: EPSG 4326 declaration (Spherical Mercator)
           r ="bilinear", # recommended method to resample continuous rasters, such as ours
           output_Raster=TRUE,
           overwrite=TRUE, 
           verbose=TRUE)
}

