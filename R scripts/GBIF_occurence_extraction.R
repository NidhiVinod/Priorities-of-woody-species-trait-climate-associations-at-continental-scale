### Script to download species' occurrences from GBIF ###
### Author: Camila Medeiros ###
### Date 12 Feb 2024 ###


##### Prep the list of species before requesting data from GBIF #### 
# Load rgbif package:
library(rgbif)

# set WD:
setwd("~/ClimateSpecies/AllSpecies")

# Load list of unique species names:
spList <- read.csv("allSpTNRS.csv", header = T)

# function to perform the name check:
name_list <- spList$Species_NameCheckedTNRS

out <- name_backbone_checklist(name_list)
write.csv(out, "SpeciesGBIFNameCheck.csv") 


##### 1. Climate- Extract species-level climate data #####
# refs: https://docs.ropensci.org/rgbif/articles/taxonomic_names.html

####  Step 1: Download species occurrences ####
# Clean R environment:
rm(list=ls())

# Load packages
library(rgbif) # to connect with GBIF and create request to download species occ data
library(dplyr) # to organize data
library(readr) # to save separate csv files for each species

# Set WD for species data:
setwd("~/ClimateSpecies/AllSpecies")

# fill in your gbif.org credentials 
user <- "" # your gbif.org username 
pwd <- "" # your gbif.org password
email <- "" # your email 

# match the accepted species names from step 0 with the gbif sp identifiers:
gbif_taxon_keys <- 
  data <- read.csv("allSpGBIF.csv") %>% #load list of accepted sp names
  #head(10) %>% # use only first 10 names for testing
  pull("Species_gbifName") %>% # use fewer names if you want to just test 
  name_backbone_checklist()  %>% # match to backbone
  filter(!matchType == "NONE") %>% # get matched names
  pull(usageKey) # get the gbif taxon keys

# create link to download all occurrences under a single DOI (!!very important here to use pred_in!!):
occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN')),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  pred_gte("year", 1950),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email)

# check download status
occ_download_wait('0062742-240506114902167')

# download occ dataframe when ready
d <- occ_download_get('0062742-240506114902167') %>%
  occ_download_import()

# save occ dataframe as .csv for future use:
write.csv(d, "NEONSpeciesOccGBIF.csv")

# Bring dataset back after download
occData <- read.csv("NEONSpeciesOccGBIF.csv")

# Check if the download includes all (and only) the species in the query:
nn <- data.frame(occData %>%  group_by(species) %>% summarise(sum = n())) 
colnames(nn) <- c("Species_gbifName", "nn") 

spList <- read.csv("allSpGBIF.csv", header = T) 

merged <- left_join(spList, nn) 


# Two of the species were sampled in 2 sites, with individuals ID'ed to the species level in one and to subspecies/variety in the other
# Thus, if using species only to subset, the subspecies disappear.
# To avoid that, here I'm creating a separate species ID for these 2 subspecies so they are not removed
Corylus_cornuta_californica <- subset(occData, taxonKey == "2876033") #using taxonKey obtaied during name check in step 0 to retrieve the subspecies
Quercus_wislizeni_frutescens <- subset(occData, taxonKey == "2878385")

Corylus_cornuta_californica$species <- "Corylus cornuta californica" #rename so they will be distinct from Corylus cornuta
Quercus_wislizeni_frutescens$species <- "Quercus wislizeni frutescens"#rename so they will be distinct from Quercus wislizeni

occDataAll <- rbind(occData, Corylus_cornuta_californica, Quercus_wislizeni_frutescens) #incorporate the subspecies back into the combined occ dataframe

# Check if we have a complete match now:
nn <- data.frame(occDataAll %>%  group_by(species) %>% summarise(sum = n())) #summarize occ by species
nn

# Save this version with all species combined:
write.csv(occDataAll, "NEONSpeciesOccGBIF.csv") 
occDataAll$group <- occDataAll$species

occDataAll %>% 
  group_by(group) %>% 
  group_walk(~ write_csv(.x, paste0("occ_raw/", .y$group, ".csv")))



####  Step 2: Clean tabular occurrence data #### 
# refs: https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_GBIF_data_with_CoordinateCleaner.html
# https://rdrr.io/cran/CoordinateCleaner/man/cc_val.html

# Clean R environment:
# rm(list=ls())

# Load packages:
library(dplyr)
library(CoordinateCleaner) 
#library(sf)
library(terra)

# Set WD for species data:
setwd("~/ClimateSpecies/AllSpecies/")

# Load species list:
spList <- read.csv("allSpGBIF.csv")
nsp <- length(spList)

for(ii in 1:nsp){ 
  # open the raw occ dataframe
  occ_raw <- read.csv(paste0("occ_raw/", spList[ii], ".csv"), skipNul = TRUE)
  
  # remove coordinates where either the latitude or longitude is missing:
  occ_noNAs <- subset(occ_raw, (!is.na(decimalLatitude)) & (!is.na(decimalLongitude))) 
  
  # remove duplicated data based on latitude and longitude:
  dups <- duplicated(occ_noNAs[c("decimalLatitude","decimalLongitude")])
  occ_unique <- data.frame(occ_noNAs[!dups,])
  
  # remove invalid coords:
  occ_unique2 <- cc_val(x = occ_unique,
                        lon = "decimalLongitude", 
                        lat = "decimalLatitude")
  
  # clean occ dataset of flagged spatial problems (incorrect/invalid/imprecise coordinates):
  flags <- clean_coordinates(x = occ_unique2, 
                             lon = "decimalLongitude", 
                             lat = "decimalLatitude",
                             countries = "countryCode",
                             species = "species",
                             tests = c("capitals", #tests a radius around adm-0 capitals. 
                                       "centroids", #tests a radius around country centroids.
                                       "equal", #tests for equal absolute longitude and latitude.
                                       "gbif",  #tests a one-degree radius around the GBIF headquarters in Copenhagen, Denmark.
                                       "zeros")) #tests for plain zeros, equal latitude and longitude and a radius around the point 0/0.
  
  # Exclude problematic records
  occ_clean <- occ_unique2[flags$.summary,]
  
  # Only keep records from field observations deposited in herbaria: 
  occ_f <- subset(occ_clean,basisOfRecord == "PRESERVED_SPECIMEN") 
  
  # Only keep records from US and contiguous countries:
  occ_f2 <- subset(occ_f, countryCode == "CA" |  countryCode == "MX" | countryCode == "US")
  
  # Limit the temporal range of data (+- 10 years of the range of Worldclim data):
  occ_f2 <- subset(occ_f2, year>=1960)
  occ_f2 <- subset(occ_f2, year<2020)
  
  # save .csvs
  write.csv(occ_f2, paste0("occ_clean/", spList[ii],".csv")) 
}


#### Step 3: Make occs spatial and clean points ####
# Load packages:
#install_github("ropensci/CoordinateCleaner")
library(dplyr)
library(CoordinateCleaner)
library(terra)
library(raster)

# Set WD:
setwd("~/ClimateSpecies/AllSpecies/")

# Load environment layers:
clim_list <- list.files("bioclim_2021", pattern=".tif$",full.names = T)
clim <- raster::stack(clim_list)

# Load species list:
spList <- read.csv("allSpGBIF.csv")
nsp <- length(spList)

# Create a lower-resolution MAP raster to use as background of sp occ maps
lowresmap <- clim[[12]]/100

for(ii in 1:nsp){
  occ_unique <- read.csv(paste0("occ_clean/", spList[ii], ".csv"))
  coordinates(occ_unique) <- ~ decimalLongitude + decimalLatitude # Assign a coordinate reference system to spatial points
  
  myCRS1 <- CRS("+init=epsg:4326") # Define the coordinate system that will be used, WGS 84
  
  crs(occ_unique) <- myCRS1 # And add Coordinate Reference System (CRS) projection
  
  # Now we need to look for and remove erroneous points/ spatial outliers:
  conditions_occ <- extract(clim, occ_unique)
  bad_records <- is.na(conditions_occ[,1])
  conditions_occ[bad_records,]
  occ_unique <- occ_unique[!bad_records,]
  
  # Spatial filter:
  cells <- cellFromXY(clim[[1]],occ_unique)
  dups <- duplicated(cells)
  occ_final <- occ_unique[!dups,]
  
  # Export occ map:
  pdf(paste0("occ_maps/", spList[ii], "_map.pdf"), width = 10, height = 7, family = "Helvetica")
  par(pin = c(10,7), mar = c(5, 5, 3, 2))
  
  plot(lowresmap, xlim = c(-168,0), ylim = c(0,80))
  plot(occ_final, add=TRUE, pch = 19, cex = .6, col = "black")
  
  dev.off()
  
  # Export shapefile:
  raster::shapefile(occ_final,
                    filename= paste0("shp_files/", spList[ii],".shp"),
                    overwrite=TRUE)}


#### Step 4: Extract environmental data for each occurrence point  ####
## Load packages:
library("tidyverse")
library("raster")
# library("dismo")
library("rgdal")

## Set WD for environmental layers:
setwd("/Volumes/CamilaBckup/climate/")

# Worldclim vars:
clim_list <- list.files("bioclim_2021", pattern=".tif$",full.names = T) # '..' leads to the path above the folder where the .rmd file is located
clim <- raster::stack(clim_list) 

# CRU AI and PET:
ai <- raster("cgiar_2022/ai_v3_yr.tif")
pet <- raster("cgiar_2022/et0_v3_yr.tif")

# Worldclim monthly variables:
month_list <- list.files("bioclim_2021_month", pattern=".tif$",full.names = T)
month <- raster::stack(month_list)

# Monthly PET:
monthPET_list <- list.files("cgiar_2023_month/Global-ET0_v3_monthly", pattern=".tif$",full.names = T)
monthPET <- raster::stack(monthPET_list)

# Monthly AI:
monthAI_list <- list.files("cgiar_2023_month/Global-AI_v3_monthly", pattern=".tif$",full.names = T)
monthAI <- raster::stack(monthAI_list)


# Load species list:
spList <- read.csv("allSpGBIF.csv")[,1]
nsp <- length(spList)

# check if a "clim" folder exists and create it if it doesn't: 
mainDir <- "~/ClimateSpecies/AllSpecies/"
subDir <- "clim/"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# loop through and extract climate data for all species:
for(ii in 1:nsp){ 
  occ_final <- rgdal::readOGR(paste0("shp_files/", spList[ii], ".shp")) 
  
  clim1 <- raster::extract(clim, occ_final) #extract bioclim data
  clim2 <- raster::extract(ai, occ_final) #extract ai data
  clim3 <- raster::extract(pet, occ_final) #extract pet data
  clim4 <- raster::extract(month, occ_final) #extract monthly bioclim data
  clim5 <- raster::extract(monthPET, occ_final) #extract monthly pet data
  clim6 <- raster::extract(monthAI, occ_final) #extract monthly ai data
  
  occ_res <- as.data.frame(cbind(clim1, clim2, clim3, clim4, clim5, clim6))
  nn <- nrow(occ_res) # get number of observations for @ species 
  occ_res$species <- rep(spList[ii], times = nn) # add column with species code
  
  write.csv(occ_res, paste0("clim/", spList[ii], "_clim.csv")) # save extracted clim data as a .csv
}


#### Step 5: Combine all species climate data into one large csv ####
# load packages:
library(dplyr)

## Set WD:
setwd("~ClimateSpecies/AllSpecies/")

# And list of all clim variables so columns are printed with the variable names:
colnames <- as.vector(read.csv("names_clim.csv", header=F)$V1)

temp <- list.files("clim", pattern="*.csv", full.names = T)
myfiles <- lapply(temp, read.delim, sep = ",") 
myfiles2 <- lapply(myfiles, setNames, colnames) 
data <- do.call("rbind", myfiles2) 
write.csv(data[,-1], "combined.csv") 


#### Step 6: Clean, organize, adjust units of the data and calculate species means ####
library(dplyr)

#  Set WD:
setwd("~/ClimateSpecies/AllSpecies/")

# Load data:
data <- read.csv("combined.csv")[,-1] 

# Adjust AI units:
data[grepl('ai', colnames(data))] <- data[grepl('ai', colnames(data))]/10000 # convert mapped values into SI units

# Remove monthly data
firstmonthlyvar <- which(colnames(data)=="ppt_jan" )
spcolumn <- which(colnames(data)=="species" )
final <- data[,-c(firstmonthlyvar:spcolumn)]
final$species <- data$species #add species column back to the end
nvars <- ncol(final) # number of environmental variables in the dataset 
final <- final %>% arrange(species)

nn <- data.frame(final %>%  group_by(species) %>% summarise(sum = n())) 

#mean
mean_sp <- aggregate(final[,-nvars], list(final$species), mean, na.rm = T)
mean_sp <- mean_sp[,-1]

#sd
sd_sp <- aggregate(final[,-nvars], list(final$species), sd, na.rm = T) 
sd_sp <- sd_sp[,-1]

#se
se_sp <- sd_sp/ sqrt(nn$sum)

# 2.5% CI
min_sp <- aggregate(final[,-nvars], list(final$species), FUN = function(x) quantile(x, probs = 0.025, na.rm = T))
min_sp <- min_sp[,-1] 

# 97.5% CI
max_sp <- aggregate(final[,-nvars], list(final$species), FUN = function(x) quantile(x, probs = 0.975, na.rm = T))
max_sp <- max_sp[,-1] 

# range
range_sp <- max_sp - min_sp # Calculate range (or "spread") of each environmental variable for each species

# cv
CV_sp <- (sd_sp/mean_sp) * 100 # Calculate CV of each environmental variable for each species

# merge
table_m <- cbind(mean_sp, se_sp, min_sp, max_sp, range_sp, CV_sp) 
table_n <- table_m[, c(matrix(1:ncol(table_m), nrow = 6, byrow = T))] 
table_f <- cbind(nn, table_n) 

# Rename columns:
aa <- c(".1", ".2", ".3", ".4", ".5")
bb <- c("_se","_min", "_max", "_range", "_CV")

for(ii in seq_along(aa)) {
  colnames(table_f) <- gsub(aa[ii], bb[ii], colnames(table_f), fixed = TRUE)
}

write.csv(table_f, "ClimateComplete.csv")


#### Step 7: Calculate Growing Season means and descriptive stats #####
library(dplyr)

# Clean R environment:
rm(list=ls())

#  Set WD:
setwd("~ClimateSpecies/AllSpecies/")

# Load list of folder names (each folder had data for one site):
data <- read.csv("combined.csv")  # combined monthly data

data[grepl('ai', colnames(data))] <- data[grepl('ai', colnames(data))]/10000 # convert mapped values into AI units


##### Separate the monthly data into individual dataframes #####
# create a folder called "GS" to store the GS monthly data:
mainDir <- "~/ClimateSpecies/AllSpecies"
subDir <- "GS/"
ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)

# create vector with list of months:
monthList <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

for(kk in 1:length(monthList)){
  month <- data[grepl(monthList[kk], colnames(data))] 
  GSthresh <- month[,3] >= 4 & month[,1] >= (2*month[,3]) # check if the climate in the occ point satisfies the GS conditions threshold (tavg is >= 4C and the ppt is >= tavg*2)
  
  climGS <-  matrix(NA, ncol = ncol(month), nrow(month))
  
  for(ii in 1:ncol(climGS)){
    for(jj in 1:nrow(climGS)){
      if(GSthresh[jj] == TRUE){ 
        climGS[jj,ii] <- month[jj,ii] 
      }
    }
  }
  climGS.df <- as.data.frame(climGS) 
  colnames(climGS.df) <- colnames(month) 
  write.csv(climGS.df, paste0("GS/", monthList[kk], ".csv")) 
}



##### Calculate the length of the potential growing season  (LPGS) #####
monthlyClim  <- list.files("GS", pattern=".csv",full.names = T) 
ppt <- matrix(NA, ncol = length(monthList), nrow = nrow(climGS)) 

for(ii in 1:length(monthList)){
  clim <- (read.csv(monthlyClim[ii]))[,2]
  ppt[,ii] <- clim 
  pptGS.df <- as.data.frame(ppt) 
  colnames(pptGS.df) <- monthList 
}

# Count the number of months that satisfy the GS conditions:
# If the GS conditions above are not met (NAs), print 0; if they are met, print 1
pptGS.df_01 <- pptGS.df
pptGS.df_01[!is.na(pptGS.df_01)] <- 1
pptGS.df_01[is.na(pptGS.df_01)] <- 0

# calculate the sum of months that satisfy the GS conditions:
data$LPGS <- rowSums(pptGS.df_01) 


##### Calculate growing season values for each occurrence point  #####
###### Growing season ppt ######
data$GSppt <- rowSums(pptGS.df, na.rm = T)

###### Growing season tavg ######
tavg <- matrix(NA, ncol = length(monthList), nrow = nrow(climGS))

for(ii in 1:length(monthList)){
  clim <- (read.csv(monthlyClim[ii]))[,4]
  tavg[,ii] <- clim
  tavgGS.df <- as.data.frame(tavg)
  colnames(tavgGS.df) <- monthList
}
data$GStavg <- rowMeans(tavgGS.df, na.rm = T)

###### Growing season pet ######
pet <- matrix(NA, ncol = length(monthList), nrow = nrow(climGS))

for(ii in 1:length(monthList)){
  clim <- (read.csv(monthlyClim[ii]))[,9]
  pet[,ii] <- clim
  petGS.df <- as.data.frame(pet)
  colnames(petGS.df) <- monthList
}
data$GSpet <- rowSums(petGS.df, na.rm = T)

###### Growing season ai ######
ai <- matrix(NA, ncol = length(monthList), nrow = nrow(climGS))

for(ii in 1:length(monthList)){
  clim <- (read.csv(monthlyClim[ii]))[,10]
  ai[,ii] <- clim
  aiGS.df <- as.data.frame(ai)
  colnames(aiGS.df) <- monthList
}
data$GSai <- rowMeans(aiGS.df, na.rm = T)


# Subset original dataframe to keep only the GS variables:
GSdataset<- data[,c(109:ncol(data))]

# Calculate number of occ for each species:
nn <- data.frame(GSdataset %>%  group_by(species) %>% summarise(sum = n()))
mean_sp <- aggregate(GSdataset[,-1], list(GSdataset$species), mean, na.rm = T)[,-1] 
sd_sp <- aggregate(GSdataset[,-1], list(GSdataset$species), sd, na.rm = T)[,-1] 
se_sp <- sd_sp/ sqrt(nn$sum)
min_sp <- aggregate(GSdataset[,-1], list(GSdataset$species), FUN = function(x) quantile(x, probs = 0.025, na.rm = T))[,-1] 
max_sp <- aggregate(GSdataset[,-1], list(GSdataset$species), FUN = function(x) quantile(x, probs = 0.975, na.rm = T))[,-1] 
range_sp <- max_sp - min_sp
CV_sp <- (sd_sp/mean_sp) * 100 

descStatsTable <- cbind(mean_sp, se_sp, min_sp, max_sp, range_sp, CV_sp) 

table_n <- descStatsTable[, c(matrix(1:ncol(descStatsTable), nrow = 6, byrow = T))] 
table_f <- cbind(nn, table_n[,c(1:ncol(table_n))])

# Rename columns:
aa <- c(".1", ".2", ".3", ".4", ".5")
bb <- c("_se","_min", "_max", "_range", "_CV")

for(ii in seq_along(aa)) {
  colnames(table_f) <- gsub(aa[ii], bb[ii], colnames(table_f), fixed = TRUE)
}

write.csv(table_f, "ClimateGS.csv")

