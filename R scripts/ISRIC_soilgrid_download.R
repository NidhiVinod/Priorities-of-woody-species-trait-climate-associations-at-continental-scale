## tutorial will show how to access SoilGrids using Web Coverage Services and the R software
## Source: https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/wcs_from_R.md

# Clean R environment:
rm(list=ls())

# Load packages:
library(XML)
library(rgdal)
library(devtools)
devtools:::install_github("gearslaboratory/gdalUtils")
library(gdalUtils)


 
## Set WD:
setwd("/Volumes/CamilaBckup/climate/ISRICsoilgrids_2022")

setwd("~/Box Sync/2_Postdoc/NEON/")

#### nitrogen ####
# define the variables for the soil property and layer of interest:
voi <- "nitrogen" # variable of interest
depth <-  "15-30cm"
quantile <-  "mean"

voi_layer <- paste(voi,depth,quantile, sep="_") # layer of interest 


# set other variables necessary for the WCS call for all kinds of requests:
wcs_path <- paste0("https://maps.isric.org/mapserv?map=/map/",voi,".map") # Path to the WCS. See maps.isric.org
wcs_service  <- "SERVICE=WCS"
wcs_version  <- "VERSION=2.0.1" # This works for gdal >=2.3; "VERSION=1.1.1" works with gdal < 2.3.


# 1. Describe the coverage layer
# First we define the request as DescribeCoverage and we create a string for the full request using also the variables previously defined
wcs_request <- "DescribeCoverage" 

wcs <- paste(wcs_path, wcs_service, wcs_version, wcs_request, sep="&")


# Then we create a XML that can be used with the gdal utility after being saved to disk:
l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", voi_layer, parent=l1)

# Save to local disk
xml.out <- "./sg.xml"
saveXML(l1, file = xml.out)

# Finally we use gdal_translate to get the geotiff locally.
gdalUtils::gdalinfo("./sg.xml")




# Set bounding box for download of TIFF (in Homolosine projection):
bb <- c(-20037500,-6729000,20037500,8600750)
igh <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"

# paste all the components of the WCS request:
wcs <- paste(wcs_path,wcs_service,wcs_version,sep="&") # This works for gdal >= 2.3


# Then we create a XML that can be used with the gdal utility after being saved to disk:
l1 <- newXMLNode("WCS_GDAL")
l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
l1.l <- newXMLNode("CoverageName", voi_layer, parent=l1)

# Save to local disk
xml.out <- "./sg.xml"
saveXML(l1, file = xml.out)



# Download raster as GeoTIFF (Warning: it can be large!)
file.out <- './test.tif'

gdal_translate(xml.out, file.out,
               tr=c(250,250), projwin=bb,
               projwin_srs =igh, co=c("TILED=YES","COMPRESS=DEFLATE","PREDICTOR=2","BIGTIFF=YES"),
               verbose=TRUE
)
