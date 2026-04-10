##Created by Nidhi Vinod
##Priorities of woody species trait-climate associations at continental scale
##Updated on April,9, 2026
##Contains the entire analysis for the MS for the species's native range data
#------------------------------------------------------------------------------------------------- 
#load libraries

dev.off()
rm(list=ls())
    
{ 
    library(tidyverse)
    library(tidyr)
    library(dplyr)
    library(dplyr)
    library(plotrix)
    library(rsq)
    library(Hmisc)
    library(corrplot)
    library(metafor)
    library(broom)
    library("paletteer")
    library(cvequality)
    library(ggplot2)
    library(car)
    library(MASS)
    library(stats)
    library(scales)
    library(patchwork)
    library(maps)
    library(raster) #for processing some spatial data
    library("rnaturalearth")
    #install.packages("rnaturalearthdata")
    library("rnaturalearthdata" ) #for downloading shapefiles
    library(sf) #for processing shapefiles
    #library(elevatr) #for downloading elevation data
    library(magrittr) #to keep things very tidy
    #library(ggspatial) #for scale bars and arrows
    library(ggpubr) #for easy selection of symbols
    library("colourpicker") #for easy selection of colors
    library(stats)
    library(devtools)
    library(viridis)
  
  

  }

#------------------------------------------------------------------------------------------------- 

#set directory      

setwd("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis")
#------------------------------------------------------------------------------------------------- 

#Correlation Matrix####

#load species data
site_climate<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
#delete empty rows
site_climate<-site_climate[-405,]

#Adding columns to a df to get a sense of col no.
# names<-colnames(site_climate)
# names<-as.data.frame(names)

#filtering data as suggested for relevant traits
new_dat<-site_climate[,c(1:4,35,41,44,47,65:71,73:76,78:80,
                         82,83,89,94:117,119:121,124,125,131,
                         134,136,138,147:149,156,157,
                         200,201,202,204,209,215,221,287)]


#remove white mountain sites because contains mostly herbacious species
dat_clean <- new_dat %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))

#filter out site details to run in correlation matrix
data<-dat_clean[,-c(1:4)]

##Correlation matrix:

#set directly where the out puts could be written
setwd("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis")

#convert data to df
data<-as.data.frame(data)#convert to dataframe

#make values absolute
data<-abs(data)

#Source Camila's code
source("~/Desktop/Continental-Scale-trait-analysis/Correlation_Matrix_code/corr_matrix_functions_23Sep2020.R") 

### Step 1: RAW
# 1.1: run raw correlation and generate 2 output files, one with the corr coefficients and one with the p-values
corr_mat(dat=data, raw=T, log=F, rank=F, 
         file1="RAWcor.csv", file2="RAWp.csv", na.action=NA)

# 1.2: summarize the results from the 2 dataframes generated above 
# output shows the raw correlation coeff followed by significance
# *p<0.05; **p<0.01; ***p<0.001; stars or ns based on the significance of the correlation.
corData <- read.csv("RAWcor.csv")[,-1] # load the RAW corr matrix minus the first column (labels)
pData <- read.csv("RAWp.csv")[,-1] # load the RAW p-value matrix minus the first column (labels)


corr_adj(corData, pData, filename="RAWcor_sig.csv") #create daframe with results

#Now let's do the same for the log and rank data!

### Step 2: LOG
# 2.1 log-transform data: 
# columns with negative values are going to be transformed as
# log10(abs(Xi)) and columns with positive and negative values are transformed as log10(Xi+abs(min(X1...Xn)))
log.data <- matrix(ncol = ncol(data), nrow = nrow(data)) # create empty matrix

for(i in 1:ncol(data)){
  if(min(data[,i], na.rm=T) <= 0 & max(data[,i],  na.rm=T) > 0){
    log.data[,i] <- log10(data[,i] - min(data[,i], na.rm=T) + 1) 
  } 
  else {log.data[,i] <- log10(abs(data[,i])) }
}

log.data <- as.data.frame(log.data) # transform into dataframe so we can change col names
colnames(log.data) <- colnames(data)

# 2.2: run raw correlation and generate 2 output files, one with the corr coefficients and one with the p-values
corr_mat(dat=log.data, raw=F, log=T, rank=F, 
         file1="LOGcor.csv",file2="LOGp.csv", na.action=NA)

# 2.3: summarize the results from the 2 dataframes generated above
corData <- read.csv("LOGcor.csv")[,-1] # load the RAW corr matrix minus the first column (labels)
pData <- read.csv("LOGp.csv")[,-1] # load the RAW p-value matrix minus the first column (labels)
corData_1<-lapply(corData,as.numeric)
corData_1 <- data.frame(matrix(unlist(corData_1), nrow=length(corData_1), byrow=FALSE))
corr_adj(corData, pData, filename="LOGcor_sig.csv") #create daframe with results


### Step 3: RANK
# 3.1 rank transform each column: 
rank.data <- matrix(ncol = ncol(data), nrow = nrow(data))

for(i in 1:ncol(data)){
  rank.data[,i] <- rank(data[,i]) # ties.method = "average"
}

rank.data <- as.data.frame(rank.data)
colnames(rank.data) <- colnames(data)

# 3.2: run raw correlation and generate 2 output files, one with the corr coefficients and one with the p-values
corr_mat(dat=rank.data, raw=F, log=F, rank=T, 
         file1="RANKcor.csv",file2="RANKp.csv", na.action=NA)

# 3.3: summarize the results from the 2 dataframes generated above
corData <- read.csv("RANKcor.csv")[,-1] # load the RAW corr matrix minus the first column (labels)
pData <- read.csv("RANKp.csv")[,-1] # load the RAW p-value matrix minus the first column (labels)

corr_adj(corData, pData, filename="RANKcor_sig.csv") #create daframe with results

### Step 4: Create corr matrix with the p-values and correlation coefficients together
RAWcor_sig <- read.csv("RAWcor_sig.csv")[,-1]
LOGcor_sig <- read.csv("LOGcor_sig.csv")[,-1]
RANKcor_sig <- read.csv("RANKcor_sig.csv")[,-1]

finalmat(RAWcor_sig,LOGcor_sig, RANKcor_sig, filename="FINALmatrix.csv")

### Step 5: Create corr matrix with YES or NO for cell highlighting step in excel
RAWp <- read.csv("RAWp.csv")[,-1]
LOGp <- read.csv("LOGp.csv")[,-1]
RANKp <- read.csv("RANKp.csv")[,-1]

high_tab(RAWp, LOGp, RANKp, filename="HighlightTable.csv")
#------------------------------------------------------------------------------------------------- 

##Highest correlation coefficient matrix for Table_S8#####

#Load data
raw<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RAWcor.csv")
rank<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RANKcor.csv")
log<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/LOGcor.csv")

#make sure the matrix are in numeric format
#load raw, rank and log
mat_1<-raw
a<-mat_1
a<-as.matrix(sapply(a, as.numeric))
a<-abs(a)

mat_2<-rank
b<-mat_2 
b<-as.matrix(sapply(b, as.numeric))
b<-abs(b)

mat_3<-log
c<-mat_3 
c<-as.matrix(sapply(c, as.numeric))
c<-abs(c)

#goal:combine all three into one matrix by selecting the highest r out of the three
#create an empty matrix
max_matrix <- matrix(0, nrow = 67, ncol = 68)

# Loop through each element to find the maximum
for (i in 1:nrow(max_matrix)) {
  for (j in 1:ncol(max_matrix)) {
    max_matrix[i, j] <- max(a[i, j], b[i, j], c[i, j])
  }
}

#convert to df
max_matrix<-as.data.frame(max_matrix)

#provide the column names
rownames(max_matrix) <- rownames(mat_1)
colnames(max_matrix)<-colnames(mat_1)

# remove the empty row
max_matrix<-max_matrix[,-1]

#add the names of the traits
max_matrix<-cbind(mat_1$X,max_matrix)

#rename column
max_matrix<-max_matrix %>% dplyr::rename("Traits"="mat_1$X")

#write the table for Table_S8
#write.csv(max_matrix,"~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/Rmax_values_si_6.csv")


##Creating a datset with highest R values####
## removes all the non-significant correlations
##Goal: calculate and inform the traits that are not correlated
raw<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RAWcor_sig.csv")
rank<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RANKcor_sig.csv")
log<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/LOGcor_sig.csv")
rawp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RAWp.csv")
rankp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RANKp.csv")
logp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/LOGp.csv")

#function for converting matrix into row and column
df<-function(data){
  data_long <- `rownames<-`(setNames(reshape(data,
                                             direction = "long",
                                             v.names = "Col",
                                             idvar = "Row",
                                             times = names(data[-1]),
                                             varying = list(names(data[-1]))), 
                                     c("Row","Col","Value")),
                            NULL)
}




#transfer the names of the columns and rows, the only purpose of p df
colnames(raw)<-colnames(rawp)
colnames(log)<-colnames(logp)
colnames(rank)<-colnames(rankp)

#add column 1, trait names
raw$X<-rawp$X
log$X<-logp$X
rank$X<-rankp$X

# #run the function for every data set
raw_lon<-df(raw)# converting matrix to rows, cols from fn above
log_lon<-df(log)
rank_lon<-df(rank)

raw_max<-raw_lon
str(raw_max)

#it doesn't matter which df we use to add the high r value
raw_max$High_R <- pmax( #gets the highest R values
  raw_lon$Value,
  log_lon$Value,
  rank_lon$Value,
  na.rm = TRUE
)
colnames(raw_max)[4] <- "Srl"

#calculate the number of sig correlations
c<-raw_max[,-c(3,4)]

climate_vars <- c("GSAI", "GSPET",
                  "GSppt", "GStavg",
                  "AI_mean", "PET_mean",
                  "MAT_mean", "MAP_mean")

# Helper: not-in operator
`%notin%` <- Negate(`%in%`)

# Keep only trait × climate correlations

corr_tc <- c %>%
  filter(
    Row %notin% climate_vars,   # keep only traits as Row
    Col %in% climate_vars       # keep only climate variables as Col
  )


#Extract significance from stars

corr_tc <- corr_tc %>%
  mutate(
    sig = str_detect(High_R, "\\*")   # TRUE if *, **, or ***
  )

trait_summary <- corr_tc %>%
  group_by(Row) %>%
  summarise(
    any_sig = any(sig),
    n_sig   = sum(sig)
  ) %>%
  ungroup()

# Keep only plant traits
plant_traits <- unique(corr_tc$Row[!(corr_tc$Row %in% climate_vars)])

# Count plant traits
length(plant_traits)  #

# Numbers for manuscript

n_traits_total <- n_distinct(trait_summary$Row)
n_traits_sig   <- sum(trait_summary$any_sig)

# Traits not correlated with any climate variable
traits_not_correlated <- trait_summary %>%
  filter(!any_sig) %>%
  pull(Row)

# View
traits_not_correlated

#removing trait-trait correlations
#follow the second column till you see gsppt, cut until there

raw_max_clean<-raw_max[-c(1:3953),]

##remove climate-climate corrs
raw_max_clean<- raw_max_clean[!(raw_max_clean$Srl %in% 60:67), ]

#separate the stars,ns

raw_max_clean <- raw_max_clean %>%
  mutate(High_R = as.character(High_R)) %>%
  # separate number from stars or "ns"
  separate(High_R,
           into = c("R", "Stars"),
           sep = "(?<=\\d)(?=(\\*+|ns)$)",
           fill = "right",
           remove = FALSE) %>%
  mutate(Number = as.numeric(R))

raw_max_clean<-raw_max_clean[,-6]

raw_max_n <- raw_max_clean %>%
  mutate(y_n = case_when(
    grepl("ns$", High_R) ~ "NO",          # ends with ns-->NO
    grepl("\\*+$", High_R) ~ "YES",       # ends with one or more *-->YES
    TRUE ~ NA_character_                 # everything else
  ))


data_sorted_sp <- raw_max_n[order(raw_max_n$Col), ]#order by climate

#rename or case when
data_sorted_sp<-data_sorted_sp %>%
  mutate(S_A = case_when(Col == 'GSppt' ~ 'GSP',
                         Col == 'GStavg' ~ 'GST',
                         Col == 'GSAI' ~ 'GSAI',
                         Col == 'GSPET' ~ 'GSPET',
                         Col == 'MAT_mean' ~ 'MAT',
                         Col == 'MAP_mean' ~ 'MAP',
                         Col == 'AI_mean' ~ 'MAAI',
                         Col == 'PET_mean' ~ 'MAPET'))



#longwise highest r for meta-analysis

write.csv(data_sorted_sp,"~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/no_hawaii_sites_highest_r.csv")

#------------------------------------------------------------------------------------------------- 

##Proportion of significant correlations#####
##Uses the highlight table from correlation matrix

#load highlight table
dat<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/HighlightTable.csv")
df<-function(data){
  data_long <- `rownames<-`(setNames(reshape(data,
                                             direction = "long",
                                             v.names = "Col",
                                             idvar = "Row",
                                             times = names(data[-1]),
                                             varying = list(names(data[-1]))), 
                                     c("Row","Col","Value")),
                            NULL)
}

# #run the function for data set
dat_lon<-df(dat)# converting matrix to rows, cols from fn above


#name the fourth col
colnames(dat_lon)[4] <- "Srl"

#removing trait-trait correlations
#you do that by following the second column till you see GSPPT, crop until there

# Define the whole lists of trait and climate variable names
trait_vars <- c(
  "e_avg_um2", "s_avg_um2",
  "d_sum_sto.mm.2", "gmax_sum_mol.m.2.s.1",
  "t_avg_tri.mm.2", "SR_unitless",
  "LA_cm2", "LMA_g.m.2",
  "LTH_mm", "LD_g.cm.3",
  "LDMC_g.g.1", "SWC_g.g.1",
  "SWMA_g.m.2", "PLAdry_.",
  "PLTHdry_.", "PA_LA",
  "PA.Drymass_cm2.g.1", "PMA_g.cm.2",
  "WD_g.cm.3", "HuberValue_unitless",
  "BarkThicknessStemDiameterRatio", "Cmass_mg.g.1",
  "Narea_g.m.2", "Nmass_mg.g.1",
  "Parea_g.m.2", "Pmass_mg.g.1",
  "Kmass_mg.g.1", "Camass_mg.g.1",
  "Mgmass_mg.g.1", "Femass_mg.g.1",
  "Bmass_mg.g.1", "Mnmass_mg.g.1",
  "Namass_mg.g.1", "Znmass_mg.g.1",
  "Cumass_mg.g.1", "Momass_mg.g.1",
  "Comass_mg.g.1", "Almass_mg.g.1",
  "Asmass_mg.g.1", "Cdmass_mg.g.1",
  "Rbmass_mg.g.1", "Srmass_mg.g.1",
  "Semass_mg.g.1", "Nimass_mg.g.1",
  "Chlarea_SPAD", "N_P",
  "C_N", "SPAD_Narea",
  "D13C_permil", "TLP_.MPa",
  "K0.1_mmol.m.2.s.1.MPa.1", "psi_K0.1_80_.Mpa",
  "Jmax_area_umol.m.2.s.1", "Vcmax_area_umol.m.2.s.1",
  "gc_mod_mol.m.2.s.1", "gc_gmax",
  "gmax_Narea_mol.g.1.s.1", "HmaxTRY_m",
  "SeedMassTRY_mg", "GSAI",
  "GSPET", "GSppt",
  "GStavg", "AI_mean",
  "PET_mean", "MAT_mean",
  "MAP_mean"
)

climate_vars <- c("GSAI", "GSPET", "GSppt", "GStavg", "AI_mean", "PET_mean", "MAT_mean", "MAP_mean")

#remove climate-climate correlations and remove those numbers in every 
#repetition
dat_lon<-dat_lon[-c(1:3953),]

##all the climate variables within Srl no. 60:67 are climate-climate correlations
dat_lon_clean<- dat_lon[!(dat_lon$Srl %in% 60:67), ]
#rename or case when
dat_lon_clean<-dat_lon_clean %>%
  mutate(S_A = case_when(Col == 'GSppt' ~ 'GSP',
                         Col == 'GStavg' ~ 'GST',
                         Col == 'GSAI' ~ 'GSAI',
                         Col == 'GSPET' ~ 'GSPET',
                         Col == 'MAT_mean' ~ 'MAT',
                         Col == 'MAP_mean' ~ 'MAP',
                         Col == 'AI_mean' ~ 'MAAI',
                         Col == 'PET_mean' ~ 'MAPET'))

#calculating yes(if significant for one/three) and no(not significant) for every cliamte variable

#add YESs into a df
yes_count <- dat_lon_clean %>%
  group_by(S_A) %>%
  summarise(Yes_Count = sum(Value == "YES"))

#add NO into a df
no_count <- dat_lon_clean %>%
  group_by(S_A) %>%
  summarise(No_Count = sum(Value == "NO"))

#bind the two
yes_no<-cbind(yes_count,no_count$No_Count)

#rename so it doesn't appear as no_count$No_Count
yes_no<-yes_no %>% rename("No_Count" = "no_count$No_Count")

#count the total
yes_no$total<-yes_no$Yes_Count+yes_no$No_Count

#the ration of YES/Total is the percentage significance
yes_no$perc_sig<-yes_no$Yes_Count/yes_no$total

#make sure it's numeric
yes_no$Yes_Count<-as.numeric(yes_no$Yes_Count)
yes_no$total<-as.numeric(yes_no$total)

# Your input data
climate_counts <- yes_no

# Define pairs to compare
pairs <- list(
  c("GSP", "GST"),
  c("GSP", "GSPET"),
  c("GSP", "GSAI"),
  c("GSP", "MAAI"),
  c("GSP", "MAPET"),
  c("GSP", "MAT"),
  c("GSP", "MAP"),
  c("GST", "GSPET"),
  c("GST", "GSAI"),
  c("GST", "MAAI"),
  c("GST", "MAPET"),
  c("GST", "MAT"),
  c("GST", "MAP"),
  c("GSPET", "GSAI"),
  c("GSPET", "MAAI"),
  c("GSPET", "MAPET"),
  c("GSPET", "MAT"),
  c("GSPET", "MAP"),
  c("GSAI", "MAAI"),
  c("GSAI", "MAPET"),
  c("GSAI", "MAT"),
  c("GSAI", "MAP"),
  c("MAAI", "MAPET"),
  c("MAAI", "MAT"),
  c("MAAI", "MAP"),
  c("MAPET", "MAT"),
  c("MAPET", "MAP"),
  c("MAP", "MAT")
)

# Function to add stars for significance
sig_stars <- function(p) {
  if (p <= 0.001) {
    return("***")
  } else if (p <= 0.01) {
    return("**")
  } else if (p <= 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Run prop.test for each pair
results <- data.frame()

for (pair in pairs) {
  group1 <- pair[1]
  group2 <- pair[2]
  
  count1 <- climate_counts[climate_counts$S_A == group1, ]
  count2 <- climate_counts[climate_counts$S_A == group2, ]
  
  test <- prop.test(
    x = c(count1$Yes_Count, count2$Yes_Count),
    n = c(count1$total, count2$total),
    alternative = "greater", 
    correct = TRUE
  )
  
  star <- sig_stars(test$p.value)
  
  results <- rbind(results, data.frame(
    Comparison = paste(group1, "vs", group2),
    Proportion1 = round(count1$Yes_Count / count1$total, 3),
    Proportion2 = round(count2$Yes_Count / count2$total, 3),
    p_value = round(test$p.value, 4),
    Significance = star
  ))
}

print(results)
results <- as.data.frame(results)
write.csv(results,"~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/prop_sig_without_hawaii.csv")
#this above analysis is the right one.

#Metawin Anaylisis#####

#load the following data
raw<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RAWcor_sig.csv")
rank<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RANKcor_sig.csv")
log<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/LOGcor_sig.csv")
rawp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RAWp.csv")
rankp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/RANKp.csv")
logp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/LOGp.csv")


#function for converting matrix into row and column
df<-function(data){
  data_long <- `rownames<-`(setNames(reshape(data,
                                             direction = "long",
                                             v.names = "Col",
                                             idvar = "Row",
                                             times = names(data[-1]),
                                             varying = list(names(data[-1]))), 
                                     c("Row","Col","Value")),
                            NULL)
}

#transfer the names of the columns and rows, the only purpose of p df
colnames(raw)<-colnames(rawp)
colnames(log)<-colnames(logp)
colnames(rank)<-colnames(rankp)

#add column 1, trait names
raw$X<-rawp$X
log$X<-logp$X
rank$X<-rankp$X

#run the function for data set
raw_lon<-df(raw)# converting matrix to rows, cols from fn above
log_lon<-df(log)
rank_lon<-df(rank)

raw_max<-raw_lon

#it doesn't matter which df we use to add the high r value
raw_max$High_R <- pmax( #gets the highest R values
  raw_lon$Value,
  log_lon$Value,
  rank_lon$Value,
  na.rm = TRUE
)
colnames(raw_max)[4] <- "Srl"

#removing trait-trait correlations
#follow the second column till you see gsppt, cut until there
# Define lists of trait and climate variable names

raw_max<-raw_max[-c(1:3953),]

#remove climate-climate correlations and remove those numbers in every 
#repetition
raw_max<- raw_max[!(raw_max$Srl %in% 60:67), ]


r2_table<-raw_max

r2_table <- r2_table %>%
  mutate(High_R = as.character(High_R)) %>%
  # separate number from stars or "ns"
  separate(High_R,
           into = c("R", "Stars"),
           sep = "(?<=\\d)(?=(\\*+|ns)$)",
           fill = "right",
           remove = FALSE) %>%
  mutate(Number = as.numeric(R))

r2_table<-r2_table %>% rename("Trait"="Row", "Climate_Variable"="Col", "R_Value"="R" )
r2_table$R_Value<-as.numeric(r2_table$R_Value)
r2_table$R2 <- r2_table$R_Value^2

# remove not significant ones

r2_table <- r2_table %>%
  arrange(Trait)
write.csv(r2_table, "~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/r_table_clean_site_clim.csv")

# Compute LogRR between predictors

climate_vars <- c("GSppt", "GStavg", "GSPET", "GSAI", "AI_mean", "PET_mean", 
                  "MAT_mean", "MAP_mean")

trait_vars <- c(
  "e_avg_um2", "i_avg_unitless", "s_avg_um2", "d_sum_sto.mm.2",
  "gmax_sum_mol.m.2.s.1", "t_avg_tri.mm.2", "LA_cm2", "LMA_g.m.2",
  "LTH_mm", "LD_g.cm.3", "LDMC_g.g.1", "LL_months", "SWC_g.g.1", "SWMA_g.m.2",
  "PLAdry_.", "PLTHdry_.", "PA_LA", "PA.Drymass_cm2.g.1", "PMA_g.cm.2", "WD_g.cm.3",
  "HuberValue_unitless", "Cmass_mg.g.1", "Nmass_mg.g.1", "Pmass_mg.g.1", "Kmass_mg.g.1",
  "Camass_mg.g.1", "Mgmass_mg.g.1", "Femass_mg.g.1", "Bmass_mg.g.1", "Mnmass_mg.g.1",
  "Namass_mg.g.1", "Znmass_mg.g.1", "Cumass_mg.g.1", "Momass_mg.g.1", "Comass_mg.g.1",
  "Almass_mg.g.1", "Asmass_mg.g.1", "Cdmass_mg.g.1", "Rbmass_mg.g.1", "Srmass_mg.g.1",
  "Semass_mg.g.1", "Nimass_mg.g.1", "Chlarea_SPAD", "N_P", "C_N", "SPAD_Narea",
  "D13C_permil", "TLP_.MPa", "K0.1_mmol.m.2.s.1.MPa.1", "psi_K0.1_80_.Mpa",
  "Jmax_area_umol.m.2.s.1", "Vcmax_area_umol.m.2.s.1", "Aarea_mod_umol.m.2.s.1",
  "gc_mod_mol.m.2.s.1", "gc_gmax", "gmax_Narea_mol.g.1.s.1", "H_m", "HmaxTRY_m",
  "SeedMassTRY_mg"
)

r2_table$R_abs <- abs(r2_table$R_Value)


pairwise <- t(combn(climate_vars, 2))
logrr_data <- data.frame()

for (i in 1:nrow(pairwise)) {
  var1 <- pairwise[i, 1]
  var2 <- pairwise[i, 2]
  for (trait in trait_vars) {
    r2_1 <- r2_table$R2[r2_table$Trait == trait & r2_table$Climate_Variable == var1]
    r2_2 <- r2_table$R2[r2_table$Trait == trait & r2_table$Climate_Variable == var2]
    
    if (length(r2_1) == 0 | length(r2_2) == 0) {
      message("Missing combination: ", trait, " with ", var1, " or ", var2)
    }
  }
}

for (i in 1:nrow(pairwise)) {
  var1 <- pairwise[i, 1]
  var2 <- pairwise[i, 2]
  for (trait in trait_vars) {
    r2_1 <- r2_table$R2[r2_table$Trait == trait & r2_table$Climate_Variable == var1]
    r2_2 <- r2_table$R2[r2_table$Trait == trait & r2_table$Climate_Variable == var2]
    
    if (length(r2_1) == 1 && length(r2_2) == 1 && !is.na(r2_1) && !is.na(r2_2) && r2_1 > 0 && r2_2 > 0) {
      logrr <- log(r2_1 / r2_2)
      logrr_data <- rbind(logrr_data, data.frame(
        Trait = trait,
        Climate_Variable_1 = var1,
        Climate_Variable_2 = var2,
        Comparison = paste0(var1, "_vs_", var2),
        LogRR = logrr,
        Variance = 0.01,
        SampleSize = nrow(r2_table)  # or whatever actual N you mean here
      ))
    }
  }
}


# Meta-analysis using metafor
results <- list()

for (cmp in unique(logrr_data$Comparison)) {
  subset_data <- subset(logrr_data, Comparison == cmp)
  res <- rma(yi = LogRR, vi = Variance, data = subset_data, method = "REML")
  results[[cmp]] <- res
  cat("\n", cmp, "\n")
  print(summary(res))
}

#  Summarize which climate variable is better
summary_table <- data.frame()

for (cmp in unique(logrr_data$Comparison)) {
  subset_data <- subset(logrr_data, Comparison == cmp)
  res <- rma(yi = LogRR, vi = Variance, data = subset_data, method = "REML")
  
  estimate <- res$b[1]
  var1 <- subset_data$Climate_Variable_1[1]
  var2 <- subset_data$Climate_Variable_2[1]
  
  interpretation <- if (estimate > 0) {
    paste(var1, "is a better predictor")
  } else if (estimate < 0) {
    paste(var2, "is a better predictor")
  } else {
    "Both are equally good predictors"
  }
  
  summary_table <- rbind(summary_table, data.frame(
    Comparison = cmp,
    LogRR_Estimate = round(estimate, 3),
    Interpretation = interpretation
  ))
}

# Print final interpretation summary
print("Summary of which variable is better predictor:")
print(summary_table)
my_r_data<-summary_table
write.csv(my_r_data,"~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/meta_win_summary_high_r_data.csv")

###testing significance
summary_table_p <- data.frame()

for (cmp in unique(logrr_data$Comparison)) {
  subset_data <- subset(logrr_data, Comparison == cmp)
  res <- rma(yi = LogRR, vi = Variance, data = subset_data, method = "REML")
  
  estimate <- res$b[1]
  pval <- res$pval[1]
  ci_lb <- res$ci.lb[1]
  ci_ub <- res$ci.ub[1]
  var1 <- subset_data$Climate_Variable_1[1]
  var2 <- subset_data$Climate_Variable_2[1]
  
  interpretation <- if (pval < 0.05) {
    if (estimate > 0) {
      paste(var1, "is a significantly better predictor (p =", round(pval, 3), ")")
    } else if (estimate < 0) {
      paste(var2, "is a significantly better predictor (p =", round(pval, 3), ")")
    } else {
      "No difference (LogRR = 0)"
    }
  } else {
    "No significant difference (p ≥ 0.05)"
  }
  
  summary_table_p <- rbind(summary_table_p, data.frame(
    Comparison = cmp,
    LogRR_Estimate = round(estimate, 3),
    CI_Lower = round(ci_lb, 3),
    CI_Upper = round(ci_ub, 3),
    P_Value = round(pval, 4),
    Interpretation = interpretation
  ))
}

# Print updated summary table
print("Summary of climate variable performance with significance:")
print(summary_table_p)
summary_table_p_my_data<-summary_table_p
write.csv(summary_table_p_my_data,"~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/metawin_for_without_hawaii_plots.csv")
# ----------------------------


#separating the R values so I can use it for the table 2 in the paper

#  Calculate mean R for each climate variable pair for Table 3####
dat<-r2_table
mean_Rs <- dat %>%
  group_by(Climate_Variable) %>%   # your climate variable column names here
  summarise(Mean_R = mean(abs(R_Value))) %>%
  ungroup()

print(mean_Rs)
write.csv(mean_Rs,"mean_r_values.csv")
#------------------------------------------------------------------------------------------------- 
#Calculating mean, SE, and SD for Table 3 climate variables####

#use the longwise Highest R data from above
yes_s<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/no_hawaii_sites_highest_r.csv")

yes_s<-yes_s
#remove rows that have NO
yes_s<- yes_s %>% 
  filter(!y_n %in% "NO")

#calulate the Mean, SD, SE for the climate variables
mean_sd_se <- yes_s %>%
  group_by(S_A) %>%  # your grouping variable
  summarise(
    mean_r = mean(abs(Number), na.rm = TRUE),
    sd_r = sd(abs(Number), na.rm = TRUE),
    se_r = sd(abs(Number), na.rm = TRUE) / sqrt(sum(!is.na(Number)))
  ) %>%
  ungroup()

write.csv(mean_sd_se, "~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis//mean_se_sd_for_without_hawaii.csv")
------------------------------------------------------------------------------------------------- 
  
##MULTIREGRESSION ANALYSIS####

#### Multiple regressions: predicting traits from MAP and MAT
#### From Camila D. Medeiros 
#### 21 Oct 25 
library(MASS)
library(MuMIn) #package used to calculate AICc
library(hier.part)
library(rdacca.hp)

# Set WD:

# Load data:
#data <- read.csv("~/Desktop/Traits-and-Climate/Raw data/trait_climate_gs_dat.csv")##without hawaii data
data <- dat_clean
#data<-data[,-1] # trait data remove the srl no. col

# Function to calculate overall model p-value:
lmp <- function (fit) {
  if (class(fit) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(fit)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


# The goal here is to predict traits from MAP and MAT and run a hierarchical partitioning analyses to quantify the contribution of each to the traits.
# So let's start by doing this for a single trait, tlp:


#### TLP 
##### a) Fit linear model usin OLS: 
fit <- lm(data$TLP_.MPa ~  data$MAP_mean + data$MAT_mean, data=data)
summary(fit) # summary of the model results
confint(fit) # confidence interval
plot(fit) # data spread plots for diagnostics
qqPlot(fit, main="QQ Plot") # qq plot

# Outputs we want to save:
summary(fit)[8] # model's multiple R2
summary(fit)[9] # model's adjusted R2
lmp(fit) # model's p-value

summary(fit)[[4]][c(1,10)] #  intercept coefficient estimate and p-value 
summary(fit)[[4]][c(2,11)] #  MAP coefficient estimate and p-value 
summary(fit)[[4]][c(3,12)] #  MAT coefficient estimate and p-value 

confint(fit)[c(1,4)] # intercept confidence interval
confint(fit)[c(2,5)] # MAP confidence interval
confint(fit)[c(3,6)] # MAT confidence interval

##### b) Hierarchical partitioning (to find out the percentage of contribution from @ predictor variable): 
pred <- as.data.frame(cbind(data$MAP_mean, data$MAT_mean))
colnames(pred) <- c("MAP", "MAT")

pred <- data.frame(MAP = data$MAP_mean, MAT = data$MAT_mean)

# Fit model
mod <- lm(TLP_.MPa ~ MAP_mean + MAT_mean, data = data)

# Run relative importance analysis
library(relaimpo)
part <- calc.relimp(mod, type = "lmg", rela = TRUE)

# Extract percent contributions
map_contrib <- part$lmg["MAP_mean"] * 100
mat_contrib <- part$lmg["MAT_mean"] * 100

# Which one is higher?
higher <- ifelse(map_contrib > mat_contrib, "MAP", "MAT")

map_contrib
mat_contrib
higher
# # Fit linear model
# mod <- lm(TLP_.MPa ~ MAP_mean + MAT_mean, data = data)
# 
# # Run relative importance (hierarchical partitioning equivalent)
# library(relaimpo)
# part <- calc.relimp(mod, type = "lmg", rela = TRUE)
# 
# # Extract % contribution of MAP and MAT
# map_contrib <- part$lmg["MAP_mean"] * 100
# mat_contrib <- part$lmg["MAT_mean"] * 100
# 
# # Print
# map_contrib
# mat_contrib




# gofs <- all.regs(data$TLP_.MPa, pred, fam = "gaussian", gof = "Rsqu", print.vars = TRUE)
# 
# part <- partition(gofs, pcan = 2, var.names = names(pred))
# 
# # Outputs we want to save:
# part$I.perc[[1]][[1]] # % contribution of MAP
# part$I.perc[[1]][[2]] # % contribution of MAT


#### All traits 
##### a) raw 
# Now that we know what we need to extract from each step, let's do it for the whole dataset.


# First, let's re-shape the dataset so it's easier to work with.
# df <- cbind(data[,1:4], data[,71], data[,70], data[,c(5:63)])
# colnames(df)[5:6] <- c("MAP_mean", "MAT_mean")
# 
# 
# # We need to extract 18 variables from the analysis, so we need a data frame with ntraits rows and 18 columns to store the results (vars+trait names column):
# out <- matrix(nrow=ncol(df), ncol=19)
# out[,1] <- colnames(df)
# 
# 
# for (i in 7:ncol(df)) {  # for each column of data = y
#   # Multiple regression
#   ff <- paste(names(df)[i], "~",  "MAP_mean + MAT_mean") # define model
#   fit <- lm(as.formula(ff), data = df) # run 
#   
#   # Hierarchical partitioning
#   pred <- as.data.frame(cbind(df$MAP_mean, df$MAT_mean)) # define predictors
#   colnames(pred) <- c("MAP", "MAT") #rename cols
#   gofs <- all.regs(df[,i], pred, fam = "gaussian", gof = "Rsqu", print.vars = TRUE)
#   part <- partition(gofs, pcan = 2, var.names = names(pred))
#   
#   out[i,2] <- summary(fit)$r.squared # model's multiple r2
#   out[i,3] <- summary(fit)$adj.r.squared # model's adjusted R2
#   out[i,4] <- lmp(fit) # model's p-value
#   out[i,5:6] <- summary(fit)[[4]][c(1,10)] #  intercept coefficient estimate and p-value 
#   out[i,7:8] <- summary(fit)[[4]][c(2,11)] #  MAP coefficient estimate and p-value 
#   out[i,9:10] <- summary(fit)[[4]][c(3,12)] #  MAT coefficient estimate and p-value 
#   out[i,11:12] <- confint(fit)[c(1,4)] # intercept confidence interval
#   out[i,13:14] <- confint(fit)[c(2,5)] # MAP confidence interval
#   out[i,15:16] <- confint(fit)[c(3,6)] # MAT confidence interval
#   out[i,17] <- part$I.perc[[1]][[1]] # % contribution of MAP
#   out[i,18] <-  part$I.perc[[1]][[2]] # % contribution of MAT
#   
#   if(out[i,17] > out[i,18]){
#     out[i,19] <-  "MAP"} else{
#       out[i,19] <- "MAT"}
# }
# 
# RAW.df<-as.data.frame(out[-c(1:6),])               # change the output matrix into a data frame (to allow naming of columns and rows)
# colnames(RAW.df)<- c("Trait","multipleR2", "adjustedR2", "model_p", "intercept", "intercept_p",
#                      "MAP_slope", "MAP_p", "MAT_slope", "MAT_p", "intercept_2.5_CI", "intercept_97.5_CI",
#                      "MAP_2.5_CI", "MAP_97.5_CI","MAT_2.5_CI", "MAT_97.5_CI",
#                      "perc_contribution_MAP", "perc_contribution_MAT", "MAP or MAT contribution higher?") # name columns
# write.csv(RAW.df, "MultipleRegressionRes_raw.csv")              # write output to text files; one for r, one for p

###Nidhi's modification to the code above because it didn't run with my R version####
#### All traits ####
##### a) raw #####

# Re-shape the dataset
df <- cbind(data[,1:4], data[,71], data[,70], data[,c(5:63)])
colnames(df)[5:6] <- c("MAP_mean", "MAT_mean")

# Initialize output matrix
out <- matrix(nrow = ncol(df), ncol = 19)
out[,1] <- colnames(df)

library(relaimpo)  # replacement for hier.part

for (i in 7:ncol(df)) {  # for each trait column
  
  # Multiple regression
  ff <- paste(names(df)[i], "~ MAP_mean + MAT_mean")
  fit <- lm(as.formula(ff), data = df)
  
  # Relative importance (replacement for hierarchical partitioning)
  part <- calc.relimp(fit, type = "lmg", rela = TRUE)
  
  # Save outputs
  out[i,2] <- summary(fit)$r.squared                     # Multiple R²
  out[i,3] <- summary(fit)$adj.r.squared                 # Adjusted R²
  out[i,4] <- lmp(fit)                                   # Model p-value
  
  out[i,5:6] <- summary(fit)[[4]][c(1,10)]               # Intercept estimate + p
  out[i,7:8] <- summary(fit)[[4]][c(2,11)]               # MAP estimate + p
  out[i,9:10] <- summary(fit)[[4]][c(3,12)]              # MAT estimate + p
  
  out[i,11:12] <- confint(fit)[c(1,4)]                   # Intercept CI
  out[i,13:14] <- confint(fit)[c(2,5)]                   # MAP CI
  out[i,15:16] <- confint(fit)[c(3,6)]                   # MAT CI
  
  # % contribution (replaces part$I.perc)
  out[i,17] <- part$lmg["MAP_mean"] * 100
  out[i,18] <- part$lmg["MAT_mean"] * 100
  
  # Which variable contributes more
  out[i,19] <- ifelse(out[i,17] > out[i,18], "MAP", "MAT")
}

# Convert to dataframe and export
RAW.df <- as.data.frame(out[-c(1:6),])
colnames(RAW.df) <- c("Trait","multipleR2", "adjustedR2", "model_p", 
                      "intercept", "intercept_p",
                      "MAP_slope", "MAP_p", "MAT_slope", "MAT_p", 
                      "intercept_2.5_CI", "intercept_97.5_CI",
                      "MAP_2.5_CI", "MAP_97.5_CI",
                      "MAT_2.5_CI", "MAT_97.5_CI",
                      "perc_contribution_MAP", "perc_contribution_MAT", 
                      "MAP or MAT contribution higher?")
write.csv(RAW.df, "MultipleRegressionRes_raw_for_without_hawaii.csv", row.names = FALSE)

####Modified code by Nidhi for log####
##### b) log #####

# Log-transform data prior to analysis:
log.data <- matrix(ncol = ncol(data), nrow = nrow(data))
for(i in 5:ncol(data)){
  if(min(data[,i], na.rm = TRUE) <= 0){
    log.data[,i] <- log10(data[,i] - min(data[,i], na.rm = TRUE) + 1) 
  } else {
    log.data[,i] <- log10(data[,i])
  }
}
log.data <- as.data.frame(log.data)
colnames(log.data) <- colnames(data) 
log.data[,1:4] <- data[,1:4]

# Re-shape the dataset
df <- cbind(log.data[,1:4], log.data[,71], log.data[,70], log.data[,c(5:63)])
colnames(df)[5:6] <- c("MAP_mean", "MAT_mean")

# Initialize output matrix
out <- matrix(nrow = ncol(df), ncol = 19)
out[,1] <- colnames(df)

library(relaimpo)  # replacement for hier.part

for (i in 7:ncol(df)) {  # for each trait column
  
  # Multiple regression
  ff <- paste(names(df)[i], "~ MAP_mean + MAT_mean")
  fit <- lm(as.formula(ff), data = df)
  
  # Relative importance (replacement for hierarchical partitioning)
  part <- calc.relimp(fit, type = "lmg", rela = TRUE)
  
  # Store outputs
  out[i,2]  <- summary(fit)$r.squared                     # Multiple R²
  out[i,3]  <- summary(fit)$adj.r.squared                 # Adjusted R²
  out[i,4]  <- lmp(fit)                                   # Model p-value
  
  out[i,5:6] <- summary(fit)[[4]][c(1,10)]                # Intercept estimate + p
  out[i,7:8] <- summary(fit)[[4]][c(2,11)]                # MAP estimate + p
  out[i,9:10] <- summary(fit)[[4]][c(3,12)]               # MAT estimate + p
  
  out[i,11:12] <- confint(fit)[c(1,4)]                    # Intercept CI
  out[i,13:14] <- confint(fit)[c(2,5)]                    # MAP CI
  out[i,15:16] <- confint(fit)[c(3,6)]                    # MAT CI
  
  # % contribution (replaces part$I.perc)
  out[i,17] <- part$lmg["MAP_mean"] * 100
  out[i,18] <- part$lmg["MAT_mean"] * 100
  
  # Which variable contributes more
  out[i,19] <- ifelse(out[i,17] > out[i,18], "MAP", "MAT")
}

# Convert to dataframe and export
LOG.df <- as.data.frame(out[-c(1:6),])
colnames(LOG.df) <- c("Trait","multipleR2", "adjustedR2", "model_p", 
                      "intercept", "intercept_p",
                      "MAP_slope", "MAP_p", "MAT_slope", "MAT_p", 
                      "intercept_2.5_CI", "intercept_97.5_CI",
                      "MAP_2.5_CI", "MAP_97.5_CI",
                      "MAT_2.5_CI", "MAT_97.5_CI",
                      "perc_contribution_MAP", "perc_contribution_MAT", 
                      "MAP or MAT contribution higher?")
write.csv(LOG.df, "MultipleRegressionRes_log.csv", row.names = FALSE)

#set working directory
setwd("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis")
##now find the best fit:
### 1. Load the two CSVs
raw <- read.csv("MultipleRegressionRes_raw_for_without_hawaii.csv")
log <- read.csv("MultipleRegressionRes_log.csv")

### 2. Make sure the trait column matches
names(raw)[1] <- "Trait"
names(log)[1] <- "Trait"

### 3. Merge the two datasets by Trait
merged <- merge(raw, log, by = "Trait", suffixes = c("_RAW", "_LOG"))

### 4. Initialize summary columns
merged$best_fit <- NA
merged$adjR2_best <- NA
merged$p_best <- NA
merged$perc_contribution_MAP_best <- NA
merged$perc_contribution_MAT_best <- NA
merged$contribution_higher_best <- NA

### 5. Determine which version (RAW or LOG) fits better for each trait
for (i in 1:nrow(merged)) {
  adj_raw <- as.numeric(merged$adjustedR2_RAW[i])
  adj_log <- as.numeric(merged$adjustedR2_LOG[i])
  
  if (!is.na(adj_raw) && !is.na(adj_log)) {
    if (adj_log > adj_raw) {
      merged$best_fit[i] <- "LOG"
      merged$adjR2_best[i] <- merged$adjustedR2_LOG[i]
      merged$p_best[i] <- merged$model_p_LOG[i]
      merged$perc_contribution_MAP_best[i] <- merged$perc_contribution_MAP_LOG[i]
      merged$perc_contribution_MAT_best[i] <- merged$perc_contribution_MAT_LOG[i]
      merged$contribution_higher_best[i] <- merged$MAP.or.MAT.contribution.higher._LOG[i]
    } else {
      merged$best_fit[i] <- "RAW"
      merged$adjR2_best[i] <- merged$adjustedR2_RAW[i]
      merged$p_best[i] <- merged$model_p_RAW[i]
      merged$perc_contribution_MAP_best[i] <- merged$perc_contribution_MAP_RAW[i]
      merged$perc_contribution_MAT_best[i] <- merged$perc_contribution_MAT_RAW[i]
      merged$contribution_higher_best[i] <- merged$MAP.or.MAT.contribution.higher._RAW[i]
    }
  }
}

### 6. Calculate overall summary stats
perc_MAP_higher <- mean(merged$contribution_higher_best == "MAP", na.rm = TRUE)
perc_MAT_higher <- mean(merged$contribution_higher_best == "MAT", na.rm = TRUE)
mean_MAP_contrib <- mean(as.numeric(merged$perc_contribution_MAP_best), na.rm = TRUE)
mean_MAT_contrib <- mean(as.numeric(merged$perc_contribution_MAT_best), na.rm = TRUE)
### 6b. Calculate SE for % contributions
se_MAP <- sd(as.numeric(merged$perc_contribution_MAP_best), na.rm = TRUE) / sqrt(sum(!is.na(merged$perc_contribution_MAP_best)))
se_MAT <- sd(as.numeric(merged$perc_contribution_MAT_best), na.rm = TRUE) / sqrt(sum(!is.na(merged$perc_contribution_MAT_best)))

### 9b. Save summary stats as a table including SE
summary_df <- data.frame(
  Results = c(
    "% traits with MAT contribution higher",
    "% traits with MAP contribution higher",
    "mean percent contribution MAT",
    "mean percent contribution MAP",
    "SE percent contribution MAT",
    "SE percent contribution MAP"
  ),
  Value = c(
    perc_MAT_higher * 100,
    perc_MAP_higher * 100,
    mean_MAT_contrib,
    mean_MAP_contrib,
    se_MAT,
    se_MAP
  )
)

write.csv(summary_df, "Merged_Raw_Log_Summary_with_SE.csv", row.names = FALSE)

### 8b. Print summary
cat("Summary:\n")
cat("% traits with MAT contribution higher:", round(perc_MAT_higher*100, 2), "%\n")
cat("% traits with MAP contribution higher:", round(perc_MAP_higher*100, 2), "%\n")
cat("Mean % contribution MAT:", round(mean_MAT_contrib, 2), "\n")
cat("Mean % contribution MAP:", round(mean_MAP_contrib, 2), "\n")
cat("SE % contribution MAT:", round(se_MAT, 2), "\n")
cat("SE % contribution MAP:", round(se_MAP, 2), "\n")

### 9. Write final table
write.csv(merged, "Merged_Raw_Log_Comparison_multiple_reg.csv", row.names = FALSE)

# ### 9b. Print summary
# cat("Summary:\n")
# cat("Mean % contribution MAP:", round(mean_MAP_contrib, 2), "\n")
# cat("Mean % contribution MAT:", round(mean_MAT_contrib, 2), "\n")
# cat("% traits with MAP contribution higher:", round(perc_MAP_higher*100, 2), "%\n")
# cat("% traits with MAT contribution higher:", round(perc_MAT_higher*100, 2), "%\n")
# 
### 10. Save summary stats as a small table
summary_df <- data.frame(
  Summary = c("mean percent contribution MAP", "mean percent contribution MAT",
              "% traits with MAP contribution higher", "% traits with MAT contribution higher",
              "SE % contribution MAT:", "SE % contribution MAP:"),
  Value = c(mean_MAP_contrib, mean_MAT_contrib, perc_MAP_higher*100, perc_MAT_higher*100,se_MAT,
            se_MAP )
)
write.csv(summary_df, "Merged_Raw_Log_Summary_multiple_reg.csv", row.names = FALSE)
#------------------------------------------------------------------------------------------------- 

##Calculate the number of correlations for the highest R#####

#load data
#data_sorted_sp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/no_hawaii_sites_highest_r.csv")
hr<-data_sorted_sp

#remove all the ns values
hr_clean <- hr[!grepl("ns$", hr$High_R), ]

hr_sep <- hr_clean %>%
  mutate(High_R = as.character(High_R)) %>%
  # separate number from stars or "ns"
  separate(High_R,
           into = c("R", "Stars"),
           sep = "(?<=\\d)(?=(\\*+|ns)$)",
           fill = "right",
           remove = FALSE) %>%
  mutate(Number = as.numeric(R))

#remove unwanted cols
hr_sep<-hr_sep[,-c(3,8)]

# make into numeric
hr_sep$R<-as.numeric(hr_sep$R)

#make the values absolute
hr_sep$R<-abs(hr_sep$R)

#Number of correlations
correlation_counts <- hr_sep %>%
  group_by(Col) %>%
  summarise(n_correlations = n()) %>%
  arrange(n_correlations)

correlation_counts

#make sure you have right directory
setwd("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis")

write.csv(correlation_counts,"without_hawaii_n_correlations.csv")
write.csv(hr_sep,"without_hawaii_r_values.csv")
#------------------------------------------------------------------------------------------------- 

###Figure 1 Species Map#####

#load data
data<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
site<-read.csv("~/Desktop/Traits-and-Climate/Raw data/archive/NEON-CTFS_02May25_sites.csv")


#remove WM, and Hawaii sites site
site <- site %>% 
  filter(!Site %in% c("White Mountain Research Center", "Laupahoehoe", "Palamanui","Palamanui_becky","Laupahoehoe_becky"))

# remove the WM, Hawaii sites in species
data <- data %>% 
  filter(!Site %in% c("White Mountain Research Center", "Laupahoehoe", "Palamanui","Palamanui_becky","Laupahoehoe_becky"))

data <- data %>%
  distinct(Species, .keep_all = TRUE)

# Example data for species
species_data <- data.frame(
  species = data$Species,
  lon = data$lon,
  lat = data$lat
)
# Example data for sites
site_data <- data.frame(
  site = site$Site_code,
  lon = site$lon,
  lat = site$lat
)

na_map <- ne_countries(scale = "medium", returnclass = "sf")
na_map <- subset(na_map, name %in% c("United States of America", "Mexico", "Canada"))

##code to make the points degrees of aridity

site_data <- data.frame(
  site = site$Site_code,
  lon = site$lon,
  lat = site$lat,
  AI = site$ai  # Ensure this column exists in your site data
)

#Use scales to have aridity legend and plot based on of AI####

# Unique sites + AI
site_ai <- unique(site_data[, c("site", "AI")])

# Build the SAME palette as the fill scale
pal <- viridis::viridis(
  256,
  option = "inferno",
  direction = -1
)

# Rescale AI to palette index
ai_rng <- range(site_ai$AI, na.rm = TRUE)
site_ai$col <- pal[
  round(rescale(site_ai$AI, to = c(1, 256), from = ai_rng))
]

# Order sites by AI 
site_ai <- site_ai[order(site_ai$AI), ]

site_levels <- site_ai$site
site_colors <- setNames(site_ai$col, site_levels)

ggplot(na_map) +
  geom_sf(fill = "lightgrey", color = "gray50", size = 0.3) +
  
  # Species points
  geom_point(
    data = species_data,
    aes(lon, lat),
    color = "steelblue", size = 3, alpha = 0.4
  ) +
  
  # Site points (fill = AI → aridity legend)
  geom_point(
    data = site_data,
    aes(lon, lat, fill = AI),
    shape = 21, color = "black", stroke = 0.7, size = 4.5
  ) +
  
  # Dummy points ONLY to create site legend
  geom_point(
    data = site_ai,
    aes(x = -999, y = -999, color = site),
    size = 4
  ) +
  
  # Aridity scale (continuous)
  scale_fill_viridis_c(
    option = "inferno",
    direction = -1,
    name = "Aridity Index"
  ) +
  
  # Site legend using colors sampled from aridity scale
  scale_color_manual(
    name = "Site",
    values = site_colors,
    breaks = site_levels
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 6),
      keyheight = unit(0.35, "cm")
    )
  )  + 
  
  coord_sf(xlim = c(-150, -40), ylim = c(22, 60), expand = FALSE) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 20),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width  = unit(0.45, "cm"),
    legend.spacing.y  = unit(0.1, "cm"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )


#------------------------------------------------------------------------------------------------- 

####Figure 2 trait-trait correlations#####
rm(list=ls())
{
  #load species data
  #dat <- read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
  dat<-dat_clean
  
  #remove y axis for individual plots except for the first one
  plot_vars <- list(
    list(xvar = "GStavg", yvar = "WD_g.cm.3", xlab = "GS T", ylab = "WD \n(g/cm³)"),
    list(xvar = "GSppt", yvar = "WD_g.cm.3", xlab = "GS PPT", ylab = ""),
    list(xvar = "GSPET", yvar = "WD_g.cm.3", xlab = "GS PET", ylab = ""),
    list(xvar = "GSAI", yvar = "WD_g.cm.3", xlab = "GS AI", ylab = ""),
    list(xvar = "MAT_mean", yvar = "WD_g.cm.3", xlab = "MAT", ylab = ""),
    list(xvar = "MAP_mean", yvar = "WD_g.cm.3", xlab = "MAP", ylab = ""),
    list(xvar = "PET_mean", yvar = "WD_g.cm.3", xlab = "PET", ylab = ""),
    list(xvar = "AI_mean", yvar = "WD_g.cm.3", xlab = "AI", ylab = "")
    
  )
  
  make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
    df <- data %>% dplyr::select(all_of(c(xvar, yvar))) %>% drop_na()
    colnames(df) <- c("x", "y")
    
    fit <- lm(y ~ x, data = df)
    summ <- summary(fit)
    r2 <- summ$r.squared
    slope <- coef(fit)[2]
    R_val <- sqrt(r2) * sign(slope)
    p_val <- summ$coefficients[2, 4]
    
    if (p_val < 0.001) {
      stars <- "***"
    } else if (p_val < 0.01) {
      stars <- "**"
    } else if (p_val < 0.05) {
      stars <- "*"
    } else {
      stars <- " NS"
    }
    
    annotation_text <- paste0("italic(r) == ", round(R_val, 2), ' * "', stars, '"')
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        plot.margin = margin(5,5,5,5),
        axis.title.y = if (ylab == "") element_blank() else element_text(),
        axis.title.x = if (xlab == "") element_blank() else element_text(),
        axis.text.x = if (xlab == "") element_blank() else element_text(),
        axis.ticks.x = if (xlab == "") element_blank() else element_line(),
        axis.text.y = if (ylab == "") element_blank() else element_text(),
        axis.ticks.y = if (ylab == "") element_blank() else element_line()
      ) +
      annotate("text",
               x = Inf, y = Inf,
               label = annotation_text,
               parse = TRUE,
               hjust = 1.1, vjust = 1.1,
               size = 3,
               fontface = "plain",   # keeps only "r" italic, not the whole thing
               color = "black")
    # Add regression line only if p-value < 0.05
    if (p_val < 0.05) {
      p <- p + geom_smooth(method = "lm", se = FALSE, color = "black")
    }
    
    return(p)
  }
  
  plots <- purrr::map(plot_vars, ~ make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat))
  
  for (i in seq_along(plots)) {
    col_num <- ((i - 1) %% 8) + 1   # 8 columns
    
    # Remove y axis elements except first column
    if (col_num != 1) {
      plots[[i]] <- plots[[i]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    # # Remove individual x-axis titles (we'll add a common one)
    # plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    
    # Keep individual x-axis titles and make x-axis text smaller
    plots[[i]] <- plots[[i]] +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text(size = 6),  # smaller x-axis tick labels
            axis.ticks.x = element_line())
  }
  final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") & 
    theme(plot.margin = margin(5,5,40,5))  # extra bottom margin for common x-axis label
  wd<-final_plot
  print(wd)
}

#WD WITHOUT X####
make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
  
  df <- data %>% 
    dplyr::select(all_of(c(xvar, yvar))) %>% 
    drop_na()
  
  colnames(df) <- c("x", "y")
  
  fit <- lm(y ~ x, data = df)
  summ <- summary(fit)
  p_val <- summ$coefficients[2, 4]
  
  stars <- if (p_val < 0.001) "***" 
  else if (p_val < 0.01) "**" 
  else if (p_val < 0.05) "*" 
  else "NS"
  
  annotation_text <- paste0(
    "italic(r) == ",
    round(sqrt(summ$r.squared) * sign(coef(fit)[2]), 2),
    ' * "', stars, '"'
  )
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(
      method = ifelse(p_val < 0.05, "lm", NA),
      se = FALSE,
      color = "black"
    ) +
    labs(x = NULL, y = ylab) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = annotation_text,
      parse = TRUE,
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  return(p)
}
plots <- purrr::map(plot_vars, ~ 
                      make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat)
)

plots[-1] <- purrr::map(
  plots[-1],
  ~ .x + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
)

final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))

wd_x <- final_plot
print(wd_x)

##D13C#####
{
  plot_vars <- list(
    list(xvar = "GStavg", yvar = "D13C_permil", xlab = "GS T", ylab = "Δ¹³C\n (‰)"),
    list(xvar = "GSppt", yvar = "D13C_permil", xlab = "GS PPT", ylab = ""),
    list(xvar = "GSPET", yvar = "D13C_permil", xlab = "GS PET", ylab = ""),
    list(xvar = "GSAI", yvar = "D13C_permil", xlab = "GS AI", ylab = ""),
    list(xvar = "MAT_mean", yvar = "D13C_permil", xlab = "MAT", ylab = ""),
    list(xvar = "MAP_mean", yvar = "D13C_permil", xlab = "MAP", ylab = ""),
    list(xvar = "PET_mean", yvar = "D13C_permil", xlab = "PET", ylab = ""),
    list(xvar = "AI_mean", yvar = "D13C_permil", xlab = "AI", ylab = "")
    
  )
  
  make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
    df <- data %>% dplyr::select(all_of(c(xvar, yvar))) %>% drop_na()
    colnames(df) <- c("x", "y")
    
    fit <- lm(y ~ x, data = df)
    summ <- summary(fit)
    r2 <- summ$r.squared
    slope <- coef(fit)[2]
    R_val <- sqrt(r2) * sign(slope)
    p_val <- summ$coefficients[2, 4]
    
    if (p_val < 0.001) {
      stars <- "***"
    } else if (p_val < 0.01) {
      stars <- "**"
    } else if (p_val < 0.05) {
      stars <- "*"
    } else {
      stars <- " NS"
    }
    
    annotation_text <- paste0("italic(r) == ", round(R_val, 2), ' * "', stars, '"')
    
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        plot.margin = margin(5,5,5,5),
        axis.title.y = if (ylab == "") element_blank() else element_text(),
        axis.title.x = if (xlab == "") element_blank() else element_text(),
        axis.text.x = if (xlab == "") element_blank() else element_text(),
        axis.ticks.x = if (xlab == "") element_blank() else element_line(),
        axis.text.y = if (ylab == "") element_blank() else element_text(),
        axis.ticks.y = if (ylab == "") element_blank() else element_line()
      ) +
      # annotate("text", x = Inf, y = -Inf, label = annotation_text, 
      #          hjust = 1.1, vjust = -1.2, size = 3, fontface = "italic", color = 'black')
      annotate("text",
               x = Inf, y = -Inf,
               label = annotation_text,
               parse = TRUE,
               hjust = 1.1, vjust = -1.2,
               size = 3,
               fontface = "plain",   # keeps only "r" italic, not the whole thing
               color = "black")
    # Add regression line only if p-value < 0.05
    if (p_val < 0.05) {
      p <- p + geom_smooth(method = "lm", se = FALSE, color = "black")
    }
    
    return(p)
  }
  
  plots <- purrr::map(plot_vars, ~ make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat))
  
  for (i in seq_along(plots)) {
    col_num <- ((i - 1) %% 8) + 1   # 8 columns
    
    # Remove y axis elements except first column
    if (col_num != 1) {
      plots[[i]] <- plots[[i]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    # # Remove individual x-axis titles (we'll add a common one)
    # plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    
    # Keep individual x-axis titles and make x-axis text smaller
    plots[[i]] <- plots[[i]] +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text(size = 6),  # smaller x-axis tick labels
            axis.ticks.x = element_line())
  }
  final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") & 
    theme(plot.margin = margin(5,5,40,5))  # extra bottom margin for common x-axis label
  
  # final_plot <- final_plot + plot_annotation(
  #   caption = "Common X Axis Label Here",
  #   theme = theme(
  #     plot.caption = element_text(size = 14, hjust = 0.5, margin = margin(t = 15))
  #   )
  # )
  d13c<-final_plot
  print(d13c)
}

#D13C WITHOUT X####

make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
  
  df <- data %>% 
    dplyr::select(all_of(c(xvar, yvar))) %>% 
    drop_na()
  
  colnames(df) <- c("x", "y")
  
  fit <- lm(y ~ x, data = df)
  summ <- summary(fit)
  p_val <- summ$coefficients[2, 4]
  
  stars <- if (p_val < 0.001) "***" 
  else if (p_val < 0.01) "**" 
  else if (p_val < 0.05) "*" 
  else "NS"
  
  annotation_text <- paste0(
    "italic(r) == ",
    round(sqrt(summ$r.squared) * sign(coef(fit)[2]), 2),
    ' * "', stars, '"'
  )
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(
      method = ifelse(p_val < 0.05, "lm", NA),
      se = FALSE,
      color = "black"
    ) +ylim(5,30)+
    labs(x = NULL, y = ylab) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = annotation_text,
      parse = TRUE,
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  return(p)
}
plots <- purrr::map(plot_vars, ~ 
                      make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat)
)

plots[-1] <- purrr::map(
  plots[-1],
  ~ .x + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
)

final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))


d13c_x<-final_plot
print(d13c_x)


##LA####
{
  plot_vars <- list(
    list(xvar = "GStavg", yvar = "LA_cm2", xlab = "GS T", ylab = "LA\n (cm²)"),
    list(xvar = "GSppt", yvar = "LA_cm2", xlab = "GS PPT", ylab = ""),
    list(xvar = "GSPET", yvar = "LA_cm2", xlab = "GS PET", ylab = ""),
    list(xvar = "GSAI", yvar = "LA_cm2", xlab = "GS AI", ylab = ""),
    list(xvar = "MAT_mean", yvar = "LA_cm2", xlab = "MAT", ylab = ""),
    list(xvar = "MAP_mean", yvar = "LA_cm2", xlab = "MAP", ylab = ""),
    list(xvar = "PET_mean", yvar = "LA_cm2", xlab = "PET", ylab = ""),
    list(xvar = "AI_mean", yvar = "LA_cm2", xlab = "AI", ylab = "")
    
  )
  
  make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
    df <- data %>% dplyr::select(all_of(c(xvar, yvar))) %>% drop_na()
    colnames(df) <- c("x", "y")
    
    fit <- lm(y ~ x, data = df)
    summ <- summary(fit)
    r2 <- summ$r.squared
    slope <- coef(fit)[2]
    R_val <- sqrt(r2) * sign(slope)
    p_val <- summ$coefficients[2, 4]
    
    if (p_val < 0.001) {
      stars <- "***"
    } else if (p_val < 0.01) {
      stars <- "**"
    } else if (p_val < 0.05) {
      stars <- "*"
    } else {
      stars <- " NS"
    }
    
    annotation_text <- paste0("italic(r) == ", round(R_val, 2), ' * "', stars, '"')
    
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      labs(x = xlab, y = ylab) +
      theme_bw() +ylim(0,300)+
      theme(
        aspect.ratio = 1,
        plot.margin = margin(5,5,5,5),
        axis.title.y = if (ylab == "") element_blank() else element_text(),
        axis.title.x = if (xlab == "") element_blank() else element_text(),
        axis.text.x = if (xlab == "") element_blank() else element_text(),
        axis.ticks.x = if (xlab == "") element_blank() else element_line(),
        axis.text.y = if (ylab == "") element_blank() else element_text(),
        axis.ticks.y = if (ylab == "") element_blank() else element_line()
      ) +
      annotate("text",
               x = Inf, y = Inf,
               label = annotation_text,
               parse = TRUE,
               hjust = 1.1, vjust = 1.1,
               size = 3,
               fontface = "plain",   # keeps only "r" italic, not the whole thing
               color = "black")
    # Add regression line only if p-value < 0.05
    if (p_val < 0.05) {
      p <- p + geom_smooth(method = "lm", se = FALSE, color = "black")
    }
    
    return(p)
  }
  
  plots <- purrr::map(plot_vars, ~ make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat))
  
  for (i in seq_along(plots)) {
    col_num <- ((i - 1) %% 8) + 1   # 8 columns
    
    # Remove y axis elements except first column
    if (col_num != 1) {
      plots[[i]] <- plots[[i]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    # # Remove individual x-axis titles (we'll add a common one)
    # plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    
    # Keep individual x-axis titles and make x-axis text smaller
    plots[[i]] <- plots[[i]] +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text(size = 6),  # smaller x-axis tick labels
            axis.ticks.x = element_line())
  }
  final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") & 
    theme(plot.margin = margin(5,5,40,5))  # extra bottom margin for common x-axis label
  
  # final_plot <- final_plot + plot_annotation(
  #   caption = "Common X Axis Label Here",
  #   theme = theme(
  #     plot.caption = element_text(size = 14, hjust = 0.5, margin = margin(t = 15))
  #   )
  # )
  la<-final_plot
  print(la)
}

#LA WITHOUT X:####
make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
  
  df <- data %>% 
    dplyr::select(all_of(c(xvar, yvar))) %>% 
    drop_na()
  
  colnames(df) <- c("x", "y")
  
  fit <- lm(y ~ x, data = df)
  summ <- summary(fit)
  p_val <- summ$coefficients[2, 4]
  
  stars <- if (p_val < 0.001) "***" 
  else if (p_val < 0.01) "**" 
  else if (p_val < 0.05) "*" 
  else "NS"
  
  annotation_text <- paste0(
    "italic(r) == ",
    round(sqrt(summ$r.squared) * sign(coef(fit)[2]), 2),
    ' * "', stars, '"'
  )
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(
      method = ifelse(p_val < 0.05, "lm", NA),
      se = FALSE,
      color = "black"
    ) +
    labs(x = NULL, y = ylab) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = annotation_text,
      parse = TRUE,
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +ylim(0,300)+
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  return(p)
}
plots <- purrr::map(plot_vars, ~ 
                      make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat)
)

plots[-1] <- purrr::map(
  plots[-1],
  ~ .x + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
)

final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))

la_x<-final_plot
print(la_x)


#NMASS####
{
  plot_vars <- list(
    list(xvar = "GStavg", yvar = "Nmass_mg.g.1", xlab = "GS T", ylab= expression(atop(italic(N[mass]), "(mg/g)"))),
    list(xvar = "GSppt", yvar = "Nmass_mg.g.1", xlab = "GS PPT", ylab = ""),
    list(xvar = "GSPET", yvar = "Nmass_mg.g.1", xlab = "GS PET", ylab = ""),
    list(xvar = "GSAI", yvar = "Nmass_mg.g.1", xlab = "GS AI", ylab = ""),
    list(xvar = "MAT_mean", yvar = "Nmass_mg.g.1", xlab = "MAT", ylab = ""),
    list(xvar = "MAP_mean", yvar = "Nmass_mg.g.1", xlab = "MAP", ylab = ""),
    list(xvar = "PET_mean", yvar = "Nmass_mg.g.1", xlab = "PET", ylab = ""),
    list(xvar = "AI_mean", yvar = "Nmass_mg.g.1", xlab = "AI", ylab = ""))
  
  make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
    df <- data %>% dplyr::select(all_of(c(xvar, yvar))) %>% drop_na()
    colnames(df) <- c("x", "y")
    
    fit <- lm(y ~ x, data = df)
    summ <- summary(fit)
    r2 <- summ$r.squared
    slope <- coef(fit)[2]
    R_val <- sqrt(r2) * sign(slope)
    p_val <- summ$coefficients[2, 4]
    
    if (p_val < 0.001) {
      stars <- "***"
    } else if (p_val < 0.01) {
      stars <- "**"
    } else if (p_val < 0.05) {
      stars <- "*"
    } else {
      stars <- " NS"
    }
    
    annotation_text <- paste0("italic(r) == ", round(R_val, 2), ' * "', stars, '"')
    
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      labs(x = xlab, y = ylab) +
      theme_bw() +ylim(0,60)+
      theme(
        aspect.ratio = 1,
        plot.margin = margin(5,5,5,5),
        axis.title.y = if (ylab == "") element_blank() else element_text(),
        axis.title.x = if (xlab == "") element_blank() else element_text(),
        axis.text.x = if (xlab == "") element_blank() else element_text(),
        axis.ticks.x = if (xlab == "") element_blank() else element_line(),
        axis.text.y = if (ylab == "") element_blank() else element_text(),
        axis.ticks.y = if (ylab == "") element_blank() else element_line()
      ) +
      annotate("text",
               x = Inf, y = Inf,
               label = annotation_text,
               parse = TRUE,
               hjust = 1.1, vjust = 1.1,
               size = 3,
               fontface = "plain",   # keeps only "r" italic, not the whole thing
               color = "black")
    # Add regression line only if p-value < 0.05
    if (p_val < 0.05) {
      p <- p + geom_smooth(method = "lm", se = FALSE, color = "black")
    }
    
    return(p)
  }
  
  plots <- purrr::map(plot_vars, ~ make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat))
  
  for (i in seq_along(plots)) {
    col_num <- ((i - 1) %% 8) + 1   # 8 columns
    
    # Remove y axis elements except first column
    if (col_num != 1) {
      plots[[i]] <- plots[[i]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    # # Remove individual x-axis titles (we'll add a common one)
    # plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    
    # Keep individual x-axis titles and make x-axis text smaller
    plots[[i]] <- plots[[i]] +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text(size = 6),  # smaller x-axis tick labels
            axis.ticks.x = element_line())
  }
  final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") & 
    theme(plot.margin = margin(5,5,40,5))  # extra bottom margin for common x-axis label
  

  nmass<-final_plot
  print(nmass)
}

#NMASS WITHOUT X####
make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
  
  df <- data %>% 
    dplyr::select(all_of(c(xvar, yvar))) %>% 
    drop_na()
  
  colnames(df) <- c("x", "y")
  
  fit <- lm(y ~ x, data = df)
  summ <- summary(fit)
  p_val <- summ$coefficients[2, 4]
  
  stars <- if (p_val < 0.001) "***" 
  else if (p_val < 0.01) "**" 
  else if (p_val < 0.05) "*" 
  else "NS"
  
  annotation_text <- paste0(
    "italic(r) == ",
    round(sqrt(summ$r.squared) * sign(coef(fit)[2]), 2),
    ' * "', stars, '"'
  )
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(
      method = ifelse(p_val < 0.05, "lm", NA),
      se = FALSE,
      color = "black"
    ) +
    labs(x = NULL, y = ylab) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = annotation_text,
      parse = TRUE,
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  return(p)
}
plots <- purrr::map(plot_vars, ~ 
                      make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat)
)

plots[-1] <- purrr::map(
  plots[-1],
  ~ .x + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
)

final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))

nmass_x<-final_plot
print(nmass_x)


##LMA_g.m.2####
#LMA_g.m.2
{
  plot_vars <- list(
    list(xvar = "GStavg", yvar = "LMA_g.m.2", xlab = "GS T", ylab = "LMA\n (g/m²)"),
    list(xvar = "GSppt", yvar = "LMA_g.m.2", xlab = "GS PPT", ylab = ""),
    list(xvar = "GSPET", yvar = "LMA_g.m.2", xlab = "GS PET", ylab = ""),
    list(xvar = "GSAI", yvar = "LMA_g.m.2", xlab = "GS AI", ylab = ""),
    list(xvar = "MAT_mean", yvar = "LMA_g.m.2", xlab = "MAT", ylab = ""),
    list(xvar = "MAP_mean", yvar = "LMA_g.m.2", xlab = "MAP", ylab = ""),
    list(xvar = "PET_mean", yvar = "LMA_g.m.2", xlab = "PET", ylab = ""),
    list(xvar = "AI_mean", yvar = "LMA_g.m.2", xlab = "AI", ylab = "")
  )
  
  make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
    df <- data %>% dplyr::select(all_of(c(xvar, yvar))) %>% drop_na()
    colnames(df) <- c("x", "y")
    
    fit <- lm(y ~ x, data = df)
    summ <- summary(fit)
    r2 <- summ$r.squared
    slope <- coef(fit)[2]
    R_val <- sqrt(r2) * sign(slope)
    p_val <- summ$coefficients[2, 4]
    
    if (p_val < 0.001) {
      stars <- "***"
    } else if (p_val < 0.01) {
      stars <- "**"
    } else if (p_val < 0.05) {
      stars <- "*"
    } else {
      stars <- " NS"
    }
    
    annotation_text <- paste0("italic(r) == ", round(R_val, 2), ' * "', stars, '"')
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      labs(x = xlab, y = ylab) +
      theme_bw() +ylim(0,300)+
      theme(
        aspect.ratio = 1,
        plot.margin = margin(5,5,5,5),
        axis.title.y = if (ylab == "") element_blank() else element_text(),
        axis.title.x = if (xlab == "") element_blank() else element_text(),
        axis.text.x = if (xlab == "") element_blank() else element_text(),
        axis.ticks.x = if (xlab == "") element_blank() else element_line(),
        axis.text.y = if (ylab == "") element_blank() else element_text(),
        axis.ticks.y = if (ylab == "") element_blank() else element_line()
      ) +
      annotate("text",
               x = Inf, y = Inf,
               label = annotation_text,
               parse = TRUE,
               hjust = 1.1, vjust = 1.1,
               size = 3,
               fontface = "plain",   # keeps only "r" italic, not the whole thing
               color = "black")
    # Add regression line only if p-value < 0.05
    if (p_val < 0.05) {
      p <- p + geom_smooth(method = "lm", se = FALSE, color = "black")
    }
    
    return(p)
  }
  
  plots <- purrr::map(plot_vars, ~ make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat))
  
  for (i in seq_along(plots)) {
    col_num <- ((i - 1) %% 8) + 1   # 8 columns
    
    # Remove y axis elements except first column
    if (col_num != 1) {
      plots[[i]] <- plots[[i]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    # # Remove individual x-axis titles (we'll add a common one)
    # plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    
    # Keep individual x-axis titles and make x-axis text smaller
    plots[[i]] <- plots[[i]] +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text(size = 6),  # smaller x-axis tick labels
            axis.ticks.x = element_line())
  }
  final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") & 
    theme(plot.margin = margin(5,5,40,5))  # extra bottom margin for common x-axis label
  lma<-final_plot
  print(lma)
}
#LMA WITHOUT X####


make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
  
  df <- data %>% 
    dplyr::select(all_of(c(xvar, yvar))) %>% 
    drop_na()
  
  colnames(df) <- c("x", "y")
  
  fit <- lm(y ~ x, data = df)
  summ <- summary(fit)
  p_val <- summ$coefficients[2, 4]
  
  stars <- if (p_val < 0.001) "***" 
  else if (p_val < 0.01) "**" 
  else if (p_val < 0.05) "*" 
  else "NS"
  
  annotation_text <- paste0(
    "italic(r) == ",
    round(sqrt(summ$r.squared) * sign(coef(fit)[2]), 2),
    ' * "', stars, '"'
  )
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(
      method = ifelse(p_val < 0.05, "lm", NA),
      se = FALSE,
      color = "black"
    ) +
    labs(x = NULL, y = ylab) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = annotation_text,
      parse = TRUE,
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5),
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 5),
      axis.ticks.x = element_line(),
      axis.line.x  = element_line()
    )
  
  return(p)
}
plots <- purrr::map(plot_vars, ~ 
                      make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat)
)

plots[-1] <- purrr::map(
  plots[-1],
  ~ .x + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
)

final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))

lma_x<-final_plot
print(lma_x)

#TLP#####

#"TLP |MPa|"
{
  plot_vars <- list(
    list(xvar = "GStavg", yvar = "TLP_.MPa", xlab = "GS T", ylab = "TLP \n(|MPa|)"),
    list(xvar = "GSppt", yvar = "TLP_.MPa", xlab = "GS PPT", ylab = ""),
    list(xvar = "GSPET", yvar = "TLP_.MPa", xlab = "GS PET", ylab = ""),
    list(xvar = "GSAI", yvar = "TLP_.MPa", xlab = "GS AI", ylab = ""),
    list(xvar = "MAT_mean", yvar = "TLP_.MPa", xlab = "MAT", ylab = ""),
    list(xvar = "MAP_mean", yvar = "TLP_.MPa", xlab = "MAP", ylab = ""),
    list(xvar = "PET_mean", yvar = "TLP_.MPa", xlab = "PET", ylab = ""),
    list(xvar = "AI_mean", yvar = "TLP_.MPa", xlab = "AI", ylab = "")
    
  )
  
  make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
    df <- data %>% dplyr::select(all_of(c(xvar, yvar))) %>% drop_na()
    colnames(df) <- c("x", "y")
    
    fit <- lm(y ~ x, data = df)
    summ <- summary(fit)
    r2 <- summ$r.squared
    slope <- coef(fit)[2]
    R_val <- sqrt(r2) * sign(slope)
    p_val <- summ$coefficients[2, 4]
    
    if (p_val < 0.001) {
      stars <- "***"
    } else if (p_val < 0.01) {
      stars <- "**"
    } else if (p_val < 0.05) {
      stars <- "*"
    } else {
      stars <- " NS"
    }
    
    annotation_text <- paste0("italic(r) == ", round(R_val, 2), ' * "', stars, '"')
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.3, color = "steelblue") +
      labs(x = xlab, y = ylab) +
      theme_bw() +
      theme(
        aspect.ratio = 1,
        plot.margin = margin(5,5,5,5),
        axis.title.y = if (ylab == "") element_blank() else element_text(),
        axis.title.x = if (xlab == "") element_blank() else element_text(),
        axis.text.x = if (xlab == "") element_blank() else element_text(),
        axis.ticks.x = if (xlab == "") element_blank() else element_line(),
        axis.text.y = if (ylab == "") element_blank() else element_text(),
        axis.ticks.y = if (ylab == "") element_blank() else element_line()
      ) +
      annotate("text",
               x = Inf, y = Inf,
               label = annotation_text,
               parse = TRUE,
               hjust = 1.1, vjust = 1.1,
               size = 3,
               fontface = "plain",   # keeps only "r" italic, not the whole thing
               color = "black")
    # Add regression line only if p-value < 0.05
    if (p_val < 0.05) {
      p <- p + geom_smooth(method = "lm", se = FALSE, color = "black")
    }
    
    return(p)
  }
  
  plots <- purrr::map(plot_vars, ~ make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat))
  
  for (i in seq_along(plots)) {
    col_num <- ((i - 1) %% 8) + 1   # 8 columns
    
    # Remove y axis elements except first column
    if (col_num != 1) {
      plots[[i]] <- plots[[i]] +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    }
    
    # # Remove individual x-axis titles (we'll add a common one)
    # plots[[i]] <- plots[[i]] + theme(axis.title.x = element_blank())
    
    # Keep individual x-axis titles and make x-axis text smaller
    plots[[i]] <- plots[[i]] +
      theme(axis.title.x = element_text(),
            axis.text.x = element_text(size = 6),  # smaller x-axis tick labels
            axis.ticks.x = element_line())
  }
  final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") & 
    theme(plot.margin = margin(5,5,40,5))  # extra bottom margin for common x-axis label
  
  tlp<-final_plot
  print(tlp)
}

##TLP WITHOUT X AXIS:#####

make_regression_plot <- function(xvar, yvar, xlab, ylab, data) {
  
  df <- data %>% 
    dplyr::select(all_of(c(xvar, yvar))) %>% 
    drop_na()
  
  colnames(df) <- c("x", "y")
  
  fit <- lm(y ~ x, data = df)
  summ <- summary(fit)
  p_val <- summ$coefficients[2, 4]
  
  stars <- if (p_val < 0.001) "***" 
  else if (p_val < 0.01) "**" 
  else if (p_val < 0.05) "*" 
  else "NS"
  
  annotation_text <- paste0(
    "italic(r) == ",
    round(sqrt(summ$r.squared) * sign(coef(fit)[2]), 2),
    ' * "', stars, '"'
  )
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(
      method = ifelse(p_val < 0.05, "lm", NA),
      se = FALSE,
      color = "black"
    ) +
    labs(x = NULL, y = ylab) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = annotation_text,
      parse = TRUE,
      hjust = 1.1, vjust = 1.1,
      size = 3
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      plot.margin = margin(5,5,5,5),
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  
  return(p)
}
plots <- purrr::map(plot_vars, ~ 
                      make_regression_plot(.x$xvar, .x$yvar, .x$xlab, .x$ylab, dat)
)

plots[-1] <- purrr::map(
  plots[-1],
  ~ .x + theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )
)

final_plot <- wrap_plots(plots, ncol = 8, nrow = 1, guides = "collect") &
  theme(plot.margin = margin(5,5,5,5))

tlp_x<-final_plot
print(tlp_x)


###Combining plots for Fig 3#####


combined <- wd_x / tlp_x / d13c_x / la_x/ nmass_x/ lma_x + 
  plot_layout(heights = c(1, 1, 1, 1,1,1)) &   # equal height for each row
  theme(
    plot.justification = "left",   # align plots to the left
  )

# Show
print(combined)

#------------------------------------------------------------------------------------------------- 

###PCA Analysis####

#### Clean R environment:
rm(list=ls())

##### Climate PCA ##### 

#load sp. data
data<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")

#cols not needed for the analysis
data<-data[,-c(3:32)]

#remove WM and Hawaii sites
data <- data %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))

# Use only the climate variables for the PCA:
clim<-data[,c("Code","MAP_mean", "MAT_mean", "AI_mean", "PET_mean", "GSAI","GSPET",
              "GSppt", "GStavg")]

clim<-as.data.frame(clim)

clim <- clim %>%
  separate(Code, into = c("Site", "Species"), sep = "_")

#log-transform data prior to analysis:
log.data <- matrix(ncol = ncol(data), nrow = nrow(data))

for(i in 3:ncol(data)){
  if(min(data[,i], na.rm=T) <= 0 & max(data[,i],  na.rm=T) > 0){
    log.data[,i] <- log10(data[,i] - min(data[,i], na.rm=T) + 1) 
  } 
  else {log.data[,i] <- log10(abs(data[,i])) }
}

log.data <- as.data.frame(log.data)
log.data[,c(1:2)] <- data[,c(1:2)]
colnames(log.data) <- colnames(data) 

# Create clim and site objects:
data_clim <- log.data[,c(3:700)]
Site <- log.data[,1]


# Re-order sites in the order we need them to be for the legend:
site_ord <- factor(Site,levels(Site)[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)])


# Rename variabless so I can display loading labels correctly:

clim_list <-data.frame(Clim=c("MAP_mean", "MAT_mean", "AI_mean", "PET_mean", "GSAI","GSPET",
                              "GSppt", "GStavg"),
                       Text=c("MAP", "MAT", "MAAI", "MAPET", "GSAI","GSPET",
                              "GSP", "GST")) 

# pca analysis:
pca<-prcomp(~ MAP_mean+ MAT_mean+ AI_mean+ PET_mean+ GSAI+ GSPET+ GSppt+ GStavg, data=log.data, center=T, scale=T) 
pc_sum<-summary(pca)
# Create a data frame from your PCA summary
# PCA summary values
std_dev <- c(2.2053, 1.3639, 0.928, 0.48031, 0.35792, 0.2249, 0.07244)
prop_var <- c(0.6079, 0.2325, 0.1076, 0.02884, 0.01601, 0.00632, 0.00066)
cum_var <- c(0.6079, 0.8405, 0.9481, 0.97694, 0.99295, 0.99927, 0.99993)

# Combine into a wide data frame
pca_table <- data.frame(
  Measure = c("Standard Deviation", "Proportion of Variance (%)", "Cumulative Proportion (%)"),
  PC1 = c(std_dev[1], prop_var[1]*100, cum_var[1]*100),
  PC2 = c(std_dev[2], prop_var[2]*100, cum_var[2]*100),
  PC3 = c(std_dev[3], prop_var[3]*100, cum_var[3]*100),
  PC4 = c(std_dev[4], prop_var[4]*100, cum_var[4]*100),
  PC5 = c(std_dev[5], prop_var[5]*100, cum_var[5]*100),
  PC6 = c(std_dev[6], prop_var[6]*100, cum_var[6]*100),
  PC7 = c(std_dev[7], prop_var[7]*100, cum_var[7]*100)
)

# Write to csv 
#write.csv(pca_table, "PCA_summary_wide.csv", row.names = FALSE)

# broken-stick plot:
plot(pca, type = "l")

# calculate and save sp scores (adjusted to a -1:1 scale): (removes empty rows)
coords <- data.frame(pca$x[,1:8]) %>%
  mutate(PC1=rescale(PC1, to = c(-1,1)), PC2=rescale(PC2,  to = c(-1,1)))

#remove the empty row in log.data
bad_rows <- which(!complete.cases(log.data[, c("MAP_mean","MAT_mean","AI_mean",
                                               "PET_mean","GSAI","GSPET",
                                               "GSppt","GStavg")]))
#empty rows:
bad_rows

#delete them
log.data<-log.data[-c(55,327,329),]

bind the two datas
scores <- cbind(log.data[1], coords)

#Write csv
#write.csv(scores, "~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/PCA_CLIM_scores_24_2_2026.csv")

# calculate loadings:
stadev<-pca$sdev

loadings<-data.frame(pca$rotation)[1:8] %>%
  mutate(PC1=rescale(PC1, to = c(-1,1)), PC2=rescale(PC2, to = c(-1,1)), 
         variable=rownames(pca$rotation), Text=clim_list$Text)

loadings$variable<-factor(loadings$variable, levels=clim_list$Clim)

# calculate corr between variable and PC axes:
var_cor_func <- function(var.loadings, comp.sdev){
  var.loadings*comp.sdev
}

var.coord  <- t(apply(pca$rotation, 1, var_cor_func, stadev))
#write.csv(var.coord, "~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/PCA_clim.corr_PCs_24_2_26.csv")


# New dataframe to store site id + sp "coordinates":
pc_sp<-cbind(log.data[,c(1:2)], coords)

# Merge trait list with loadings:
mutate(data.frame(pca$rotation)[1:2], variable=rownames(pca$rotation), Text=clim_list$Text) %>%
  arrange(desc(abs(PC2)))

# % variance explained by each PCA axis and total (PC1+PC2):
pc1_var <-round(summary(pca)$importance[2,1]*100, 1)
pc2_var <-round(summary(pca)$importance[2,2]*100, 1)
pc3_var <-round(summary(pca)$importance[2,3]*100, 1)
pc4_var <-round(summary(pca)$importance[2,4]*100, 1)
pc5_var <-round(summary(pca)$importance[2,5]*100, 1)
pc6_var <-round(summary(pca)$importance[2,6]*100, 1)
pc7_var <-round(summary(pca)$importance[2,7]*100, 1)
pc8_var <-round(summary(pca)$importance[2,8]*100, 1)
#pc9_var <-round(summary(pca)$importance[2,9]*100, 1)

## PCA Plot:
# setup for PCA plot:
colors <- rep("black", 6)
# sites <-  c(rep("lightsteelblue1", 21), rep("cadetblue1", 9) ,rep("skyblue3",25),
#             rep("deepskyblue2", 17), rep("royalblue2", 14), rep("navyblue", 21))
# shp <- c(rep(21, 21), rep(22, 9) ,rep(23, 25),
#          rep(24, 17), rep(25, 14), rep(21, 21))


dev.off()

loadings$variable <- rownames(loadings)

##PCA (Figure 1c) Aligned based on aridity####

#make sure these libraries are loaded

#load data
site<-read.csv("~/Desktop/Traits-and-Climate/Raw data/archive/NEON-CTFS_02May25_sites.csv")

site <- site %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))

# Separate Site and Species
dat <- data %>% separate(Code, into = c("Site", "Species"), sep = "_")
dat<-site
pc_sp <- pc_sp %>% separate(Code, into = c("Site", "Species"), sep = "_")

#PCA Figure 1c####

# Get AI per site
site_ai <- dat %>%
  dplyr::select(Site_code, ai) %>%
  distinct(Site_code, .keep_all = TRUE)

# Build the SAME palette used everywhere
pal <- viridis::viridis(
  256,
  option = "inferno",
  direction = -1
)

# Rescale AI to palette indices
ai_rng <- range(site_ai$ai, na.rm = TRUE)

site_ai$color <- pal[
  round(rescale(site_ai$ai, to = c(1, 256), from = ai_rng))
]

# Order sites wet → dry (optional but recommended)
site_ai <- site_ai %>% arrange(ai)

# Named vector for ggplot
site_levels_ordered <- site_ai$Site_code
site_colors <- setNames(site_ai$color, site_levels_ordered)
pc_sp$Site <- factor(pc_sp$Site, levels = site_levels_ordered)

ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  
  # PCA points (site-colored, but aridity-derived)
  geom_point(
    data = pc_sp,
    aes(PC1, PC2, fill = Site),
    shape = 21, color = "black", size = 5
  ) +
  
  # Loadings
  geom_segment(
    data = loadings,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.03, "npc")),
    color = "black"
  ) +
  
  geom_label(
    data = loadings,
    aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Text),
    size = 3.5,
    fontface = "bold",
    color = "black",
    fill = alpha("white", 0.6),
    label.size = 0
  ) +
  
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)")
  ) +
  
  coord_equal() +
  scale_fill_manual(values = site_colors) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text = element_text(size = 19),
    axis.title = element_text(size = 19),
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )
#------------------------------------------------------------------------------------------------- 

# Latitude vs. PCA plot SI Figure 1a and b####
#Can we please test the correlation of PC1 and PC2 with latitude? 
#   I wonder if PC2 captures latitudinal variation and thus winter cold
site<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")

site <- site %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))


#use the data$lat to draw correlation from PC1 and PC2 with lat and lon
pca_sp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/PCA_CLIM_scores_24_2_2026.csv")
pca_clim<-read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/PCA_clim.corr_PCs_24_2_26.csv")

##column bind pca 1 and 2 of species with lat and long
#species data that has the lat and long. the above data doesn't have lat and long
data<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
data <- data %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))


merged_data <- left_join(pca_sp, data, by = "Code")
ggplot(data=merged_data, aes(y=PC2, x=lat))+geom_point()+geom_smooth(method = 'lm',,color="blue")+
  xlab("Latitude")+theme_bw()+ylab("Climate-PC2")

model <- lm(PC1 ~ lat, data = merged_data)
summary(model)

#edited code:
label_text <- paste0("italic(r) == ", sprintf("%.2f", r_val), " * '***'")
merged_data <- left_join(pca_sp, data, by = "Code")

model <- lm(PC1 ~ lat, data = merged_data)
r_val <- sqrt(summary(model)$r.squared)

label_text <- paste0("italic(r) == ", sprintf("%.2f", r_val), " * '***'")

ggplot(data = merged_data, aes(y = PC1, x = lat)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  xlab("Latitude") +
  ylab("Climate-PC1") +
  theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    size= 5,
    parse = TRUE,
    hjust = 1.1, vjust = 1.5
  )

##PC2
label_text <- paste0("italic(r) == ", sprintf("%.2f", r_val), " * '***'")
merged_data <- left_join(pca_sp, data, by = "Code")

model <- lm(PC2 ~ lat, data = merged_data)
r_val <- sqrt(summary(model)$r.squared)

label_text <- paste0("italic(r) == ", sprintf("%.2f", r_val), " * '***'")

ggplot(data = merged_data, aes(y = PC2, x = lat)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  xlab("Latitude") +
  ylab("Climate-PC2") +
  theme_bw() +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    size= 5,
    parse = TRUE,
    hjust = 1.1, vjust = 1.5
  )


# Residual standard error: 0.191 on 359 degrees of freedom
# Multiple R-squared:  0.431,	Adjusted R-squared:  0.4294 
# F-statistic:   272 on 1 and 359 DF,  p-value: < 2.2e-16

ggplot(data=merged_data, aes(y=PC2, x=lat))+geom_point()+geom_smooth(method = 'lm',,color="blue")+
  xlab("Latitude")+ylab("Climate-PC2")+theme_bw()

model <- lm(PC1 ~ lat, data = merged_data)## change PC1 or 2 based on the plot
summary(model)
#------------------------------------------------------------------------------------------------- 
#Figure 4#### R mean values with SE####
####color palatte used in the MS

r_dat <- read.csv("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/mean_se_sd_for_without_hawaii.csv")

r_dat<-r_dat%>%rename("Climate" = "S_A")

# Add MA vs GS grouping and CLD letters
r_dat <- r_dat %>%
  mutate(
    GS_MA = ifelse(grepl("^GS", Climate), "GS", "MA"),
    CLD = case_when(
      Climate == "MAAI"  ~ "A",
      Climate == "MAP"   ~ "A",
      Climate == "MAPET" ~ "A",
      Climate == "MAT"   ~ "B",
      Climate == "GSP"   ~ "a",
      Climate == "GST"   ~ "ab",
      Climate == "GSPET" ~ "ab",
      Climate == "GSAI"  ~ "b"
    )
  ) %>%
  mutate(
    GS_MA = factor(GS_MA, levels = c("MA", "GS"))
  ) %>%
  arrange(GS_MA, mean_r) %>%
  mutate(
    Climate = factor(Climate, levels = Climate)
  )

#Set1 + orange palette 
ma_colors <- c(
  "MAT"   = "#D73027",  # red
  "MAP"   = "#4575B4",  # blue
  "MAPET" = "#FF7F00",  # bright orange
  "MAAI"  = "#762A83"   # purple
)

gs_colors <- alpha(ma_colors, 0.75)  # lighter / semi-transparent

climate_colors <- c(ma_colors, gs_colors)
names(climate_colors) <- c("MAT","MAP","MAPET","MAAI","GST","GSP","GSPET","GSAI")

# Plot
ggplot(r_dat, aes(x = Climate, y = mean_r, fill = Climate)) +
  
  geom_bar(stat = "identity", width = 0.75) +
  
  geom_errorbar(aes(ymin = mean_r - se_r, ymax = mean_r + se_r), width = 0.15) +
  
  geom_text(aes(label = CLD, y = mean_r + se_r + 0.025),
            size = 5, fontface = "bold") +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 17),
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  ) +
  
  ylab("Mean r values") +
  ylim(0, 0.5) +
  
  scale_fill_manual(values = climate_colors, name = "Climate Group")

#------------------------------------------------------------------------------------------------- 

#ANOVA Table S_4: Site and within species####

# Load packages:
library(car)
library(MASS)
library(stats)
library(dplyr)
library(tidyverse)
## Clean R environment:
rm(list=ls())

#load data
data_tr<-read.csv("~/Desktop/Traits-and-Climate/Raw data/individuals_16Feb26.csv")
#data<-data[,c(3,4,51,52,53)]
#remove white mountain sites because it mostly has herbacious species
dat_clean <- data_tr %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))


#first chunk of code from cormat above
#load data
site_climate<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
#delete empty row
site_climate<-site_climate[-405,]

site_climate <- site_climate %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))

#filtering trait data
new_dat<-site_climate[,c(1:4,35,41,44,47,65:71,73:76,78:80,
                         82,83,89,94:117,119:121,124,125,131,
                         134,136,138,147:149,156,157,
                         200,201,202,204,209,215,221,287)]


# data_trim <- dat_clean %>%
#   dplyr::select(any_of(names(new_dat)))
new_dat<-new_dat[,-c(64:71)]#removing climate variables
data<-new_dat


## LOG transform traits:
data.l <- matrix(ncol = ncol(data), nrow = nrow(data))

for(i in 5:ncol(data)){#start the 5: where the traits start
  if(min(data[,i], na.rm=T) <= 0 & max(data[,i],  na.rm=T) > 0){
    data.l[,i] <- log10(data[,i] - min(data[,i], na.rm=T) + 1) 
  } 
  else {data.l[,i] <- log10(abs(data[,i])) }
}

data.l <- as.data.frame(data.l)
colnames(data.l) <- colnames(data) 

data.l[,c(1:4)] <- data[,c(1:4)]## 1:6 is cateogorical in camila's


sp <- data.l$Species
site <- data.l$Site
traits <- data.l[,c(5:ncol(data.l))]

df <- as.data.frame(cbind(sp,site,traits))


# aov function runs one-way ANOVA:
ntrt <- ncol(df)
anova.res <- matrix(NA, nrow= 13, ncol= ntrt) #because we are extracting 13 parameters from the anova res tables

anovares <- aov(df[, 3] ~ site + Error(sp), data = df)
print(summary(anovares))


###Sometimes this chunk has problems depending on R version
for(j in 1:13){
  for(i in 3:ntrt){
    anovares <- aov( df[,i]~ site + Error(sp), data = df) 
    anova.res[1,i] <- summary(anovares)[[1]][[1]][[1]][[1]] #df sites
    anova.res[2,i] <- summary(anovares)[[1]][[1]][[1]][[2]] #df sp
    anova.res[3,i] <- summary(anovares)[[2]][[1]][[1]][[1]] #df indv
    anova.res[4,i] <- summary(anovares)[[1]][[1]][[2]][[1]] #Sum sq sites
    anova.res[5,i] <- summary(anovares)[[1]][[1]][[2]][[2]] #Sum sq residuals
    anova.res[6,i] <- summary(anovares)[[2]][[1]][[2]][[1]] #Sum sq error
    anova.res[7,i] <- summary(anovares)[[1]][[1]][[3]][[1]] #Mean sq between sites
    anova.res[8,i] <- summary(anovares)[[1]][[1]][[3]][[2]] #Mean Sq within sites
    anova.res[9,i] <- summary(anovares)[[2]][[1]][[3]][[1]]#Mean Sq within sp
    anova.res[10,i] <- summary(anovares)[[1]][[1]][[4]][[1]] #F between sites
    anova.res[11,i] <- summary(anovares)[[1]][[1]][[5]][[1]] #p between sites
    anova.res[12,i] <- summary(anovares)[[1]][[1]][[3]][[2]]/summary(anovares)[[2]][[1]][[3]][[1]]#F within sites
    anova.res[13,i] <- pf(q= summary(anovares)[[1]][[1]][[3]][[1]]/summary(anovares)[[2]][[1]][[3]][[1]],
                          df1=summary(anovares)[[1]][[1]][[1]][[2]],
                          df2=summary(anovares)[[1]][[1]][[1]][[1]],
                          lower.tail=F)  #p within sites
  }
}


colnames(anova.res) <- colnames(df)
rownames(anova.res) <- c("df.sites", "df.sp", "df.indv", 
                         "SumSqSection", "SumSqResiduals", "SumSqError", 
                         "MeanSqbetweenSites", "MeanSqwithinSites", "MeanSqwithinSp",
                         "FbetweenSites","PbetweenSites", 
                         "FwithinSites","PwithinSites")
t<-anova.res
write.csv(anova.res, file = "ANOVA_LOG_rescomplete_jan42026.csv") #LOG


### Tukey test:
library(TukeyC)

fit <- aov(LA_cm2 ~ site + Error(sp), data=df) #need to do this one trait at a time for now
summary(fit)

#Tukey means test:
tuk <- TukeyC(fit,  data = 'df', error = 'sp', which = 'site', sig.level = 0.05)
summary(tuk)


### Table of summarized results:
df <- t(anova.res)[-c(1:2),]
restab <- matrix(NA, nrow= nrow(df), ncol= 4)

for(ii in 1:nrow(df)){
  restab[,1] <- row.names(df)
  restab[,2] <- paste0(round(df[,7],3), ", ", df[,1], ", ", round(df[,11],3), ", ", round((df[,4]/(df[,4] + df[,5]+ df[,6])*100),0), "%") #Mean sq sites, df sites, p-value, % variation explained by between sites variation
  restab[,3] <- paste0(round(df[,8],3), ", ", df[,2], ", ", round(df[,13],3), ", ", round((df[,5]/(df[,4] + df[,5]+ df[,6])*100),0), "%") #Mean sq within sites, df sp, p-value, % variation explained by between species variation
  restab[,4] <- paste0(round(df[,9],3), ", ", df[,3], ", ", round((df[,6]/(df[,4] + df[,5]+ df[,6])*100),0), "%") #Mean sq within sp, df sp, % variation explained by within sp variation
}

colnames(restab) <- c("Trait", "Site", "Species (within-site variance)", "Error (within-species variance)")

write.csv(restab,"/Users/nidhivinod/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/anova_results_feb25_2026.csv")

## ANOVA for across species and error####

rm(list=ls())

# Load data
data_tr <- read.csv("~/Desktop/Traits-and-Climate/Raw data/individuals_16Feb26.csv")
dat_clean <- data_tr %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))

site_climate <- read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
#site_climate <- site_climate[-405,]
site_climate <- site_climate %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))


new_dat<-site_climate[,c(1:4,35,41,44,47,65:71,73:76,78:80,
                         82,83,89,94:117,119:121,124,125,131,
                         134,136,138,147:149,156,157,
                         200,201,202,204,209,215,221,287)]

# data_trim <- dat_clean %>% dplyr::select(any_of(names(new_dat)))
# data <- data_trim

new_dat<-new_dat[,-c(64:71)]
data<-new_dat
# Log-transform traits (from 5th column onwards)
data.l <- matrix(ncol = ncol(data), nrow = nrow(data))

for(i in 5:ncol(data)){
  if(min(data[,i], na.rm=TRUE) <= 0 & max(data[,i], na.rm=TRUE) > 0){
    data.l[,i] <- log10(data[,i] - min(data[,i], na.rm=TRUE) + 1) 
  } else {
    data.l[,i] <- log10(abs(data[,i]))
  }
}

data.l <- as.data.frame(data.l)
colnames(data.l) <- colnames(data) 
data.l[,1:4] <- data[,1:4]

sp <- as.factor(data.l$Species)
site <- as.factor(data.l$Site)
traits <- data.l[,5:ncol(data.l)]

df <- data.frame(sp = sp, site = site, traits)

ntrt <- ncol(df)

# Create matrix to store ANOVA results
# Rows: for each metric for Site, Species(within-site), Error, Across Species (Species main effect)
# Columns: one per trait (+ 2 for sp and site labels)
# We'll store: F-value, df, p-value, % variance explained
# To store all info, create a list and then make a dataframe later

results_list <- list()

for (i in 3:ntrt) {
  
  trait_name <- colnames(df)[i]
  
  fit <- aov(
    as.formula(paste(trait_name, "~ sp")), # across species regarldless of site
    data = df
  )
  
  summ <- summary(fit)[[1]]
  
  SS_sp <- summ["sp", "Sum Sq"]
  df_sp <- summ["sp", "Df"]
  MS_sp <- summ["sp", "Mean Sq"]
  F_sp  <- summ["sp", "F value"]
  p_sp  <- summ["sp", "Pr(>F)"]
  
  SS_res <- summ["Residuals", "Sum Sq"]
  df_res <- summ["Residuals", "Df"]
  MS_res <- summ["Residuals", "Mean Sq"]
  
  total_SS <- SS_sp + SS_res
  
  var_sp <- SS_sp / total_SS * 100
  var_err <- SS_res / total_SS * 100
  
  sp_str <- paste0(
    round(F_sp, 3), ", ",
    df_sp, ", ",
    signif(p_sp, 3), ", ",
    round(var_sp, 0), "%"
  )
  
  err_str <- paste0(
    round(MS_res, 3), ", ",
    df_res, ", ",
    round(var_err, 0), "%"
  )
  
  results_list[[trait_name]] <- c(sp_str, err_str)
}

# Convert list to dataframe
results_df <- do.call(rbind, results_list)
colnames(results_df) <- c(
  "Across species",      # F, df, p, %variance
  "Residual (error)"
)  
results_df <- data.frame(Trait = rownames(results_df), results_df, row.names = NULL)

# Write out the results
write.csv(results_df, "~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis/anova_results_across_species_try_feb_26.csv", row.names = FALSE)

# View example
print(results_df)
#--------------------------------------------------------------------------------------------------------------------------
##Pie chart for Figure 3c#### 
#The numbers are from multiple regression
values <- c(MAP = 77.631822, MAT = 22.368178)

colors <- c(
  "MAP" = "#4575B4", # blue
  "MAT" = "#D73027")  # red)  

pie3D(values,
      labels = NA,          # no names or values
      explode =  0.1, # slightly separate slices
      theta = 1.5,          # faces more toward screen
      col = colors)
#------------------------------------------------------------------------------------------------- 
##CV, range, slrt and leven's test for mat and map####

#load data
site_climate<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
#delete empty row
site_climate<-site_climate[-405,]

#such a huge dataset that I am putting the colnames into a df to see which to remove
# names<-colnames(site_climate)
# names<-as.data.frame(names)

#filtering data as suggested by Lawren
new_dat<-site_climate[,c(1:4,35,41,44,47,65:71,73:76,78:80,
                         82,83,89,94:117,119:121,124,125,131,
                         134,136,138,147:149,156,157,
                         200,201,202,204,209,215,221,287)]


#remove white mountain sites because it mostly has herbacious species
dat_clean <- new_dat %>% 
  filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))

# CV function
cv <- function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)}

# Standard error of CV (Feltz & Miller 1996 formula)
cv_se <- function(x) {
  n <- sum(!is.na(x))
  cv_x <- cv(x)
  cv_x * sqrt((1 + 0.5 * cv_x^2) / n)
}

slrt_cv <- function(x, y) {
  cv1 <- cv(x)
  cv2 <- cv(y)
  se1 <- cv_se(x)
  se2 <- cv_se(y)
  
  z <- (cv1 - cv2) / sqrt(se1^2 + se2^2)
  p_value <- 2 * (1 - pnorm(abs(z)))
  
  data.frame(
    CV1 = cv1,
    CV2 = cv2,
    SE1 = se1,
    SE2 = se2,
    Z = z,
    p_value = p_value
  )
}

# raw MAP vs MAT
map <- dat_clean$MAP_mean
mat <- dat_clean$MAT_mean
map_mat_slrt<-slrt_cv(map, mat)

# log(+2) MAP vs log(+2) MAT
log_map <- log(map + 2)
log_mat <- log(mat + 2)
log_map_mat_slrt<-slrt_cv(log_map, log_mat)

map_mat_slrt<-map_mat_slrt%>% rename(Map_cv="CV1", Mat_cv="CV2",
                                     Map_se= "SE1", Mat_se="SE2" )
log_map_mat_slrt<-log_map_mat_slrt%>% rename(log_Map_cv="CV1", log_Mat_cv="CV2",
                                             log_Map_se= "SE1", log_Mat_se="SE2" )

##SLRT test:
slrt_true <- function(x, y) {
  
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  
  n1 <- length(x)
  n2 <- length(y)
  
  m1 <- mean(x)
  m2 <- mean(y)
  
  s1 <- sd(x)
  s2 <- sd(y)
  
  cv1 <- s1 / m1
  cv2 <- s2 / m2
  
  # Unconstrained log-likelihood
  ll1 <- sum(dnorm(x, m1, s1, log = TRUE)) +
    sum(dnorm(y, m2, s2, log = TRUE))
  
  # Under H0: equal CV
  # sigma = CV * mean
  # need to estimate common CV
  
  cv_pool <- (cv1*n1 + cv2*n2)/(n1+n2)
  
  s1_0 <- cv_pool * m1
  s2_0 <- cv_pool * m2
  
  ll0 <- sum(dnorm(x, m1, s1_0, log = TRUE)) +
    sum(dnorm(y, m2, s2_0, log = TRUE))
  
  LR <- 2 * (ll1 - ll0)
  
  SLRT <- sign(cv1 - cv2) * sqrt(LR)
  
  p_value <- 2 * (1 - pnorm(abs(SLRT)))
  
  data.frame(
    CV1 = cv1,
    CV2 = cv2,
    SLRT = SLRT,
    p_value = p_value
  )
}

slrt_true(map, mat)
slrt_true(log_map, log_mat) 

### Range 

library(car)   # for Levene's test

# Extract MAP and MAT
dat_clean2<-dat_clean[-c(327),]
map <- dat_clean2$MAP_mean
mat <- dat_clean2$MAT_mean

#calculate the range

range_map<-map/mean(map)
range_mat <- mat/mean(mat)

# Compute SD
sd_map <- sd(range_map, na.rm = TRUE)
sd_mat <- sd(range_mat, na.rm = TRUE)

sd(mat)
sd(map)




# Create a data frame
df <- data.frame(
  value = c(range_map, range_mat),
  group = factor(rep(c("MAP", "MAT"), c(length(range_map), length(range_mat))))
)

# Run Levene's test (center = median is robust)
levene_test <- leveneTest(value ~ group, data = df, center = median)

levene_test

##for logged data:

log_mat <- log(mat + 2)
log_map <- log(map + 2)

sd(log_map)
sd(log_mat)

#calculate the range

range_log_map<-log_map/mean(log_map)
range_log_mat <- log_mat/mean(log_mat)

# Compute SD
sd_log_map <- sd(range_log_map, na.rm = TRUE)
sd_log_mat <- sd(range_log_mat, na.rm = TRUE)


# Create a data frame
df <- data.frame(
  value = c(range_log_map, range_log_mat),
  group = factor(rep(c("MAP", "MAT"), c(length(range_log_map), length(range_log_mat))))
)

# Run Levene's test (center = median is robust)
levene_test <- leveneTest(value ~ group, data = df, center = median)

levene_test

##calculate number of species per trait####
setwd("/Users/nidhivinod/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis")
dat_clean$Species[dat_clean$Species == "Carya glabra var. odorata"] <- "Carya glabra"
dat_clean$Species[dat_clean$Species == "Cornus nuttallii"] <- "Cornus nuttalii"
dat_clean$Species[dat_clean$Species == "Fraxinus pennsylvannica"] <- "Fraxinus pennsylvanica"
dat_clean$Species[dat_clean$Species == "Prunus serrotina"] <- "Prunus serotina"
dat_clean$Species[dat_clean$Species == "Quercus wislizeni"] <- "Quercus wislizeni var. frutescens"

species_per_trait <- dat_clean %>%
  summarise(across(
    -Species,
    ~ n_distinct(Species[!is.na(.)])
  )) %>%
  pivot_longer(everything(),
               names_to = "Trait",
               values_to = "Number_of_Species")

dat_clean %>%
  filter(is.na(GStavg)) %>%
  distinct(Species)

write.csv(species_per_trait, "Number_of_sp_per_trait.csv")


## Table S_11####
# One-way ANOVA on log(R_absolute_values)

# install if needed
# install.packages(c("tidyverse", "car", "emmeans", "multcompView"))

library(tidyverse)
library(car)
library(emmeans)
library(multcompView)
library(agricolae)

#load data
df <- read.csv("without_hawaii_r_values.csv")


# log-transform r-values if needed
dat <- df %>%
  mutate(logR = log(R + 1e-6))  # +1e-6 to avoid log(0)

dat<-dat%>%rename ("ClimateVar"="Col")

# ANOVA on all GS variables together####
dat <- dat %>%
  mutate(Climate = case_when(
    str_starts(ClimateVar, "GS") ~ "GS",
    str_ends(ClimateVar, "_mean") ~ "MeanAnnual",
    #ClimateVar %in% c("MAAI") ~ "MeanAnnual",
    TRUE ~ "Other"
  ))

#log-transform r-values
dat_all <- dat %>% mutate(logR = log(R) + 1e-6)

#One-way ANOVA
anova_model <- aov(logR ~ ClimateVar, data = dat_all)

#ANOVA summary
anova_summary <- summary(anova_model)
anova_summary

#extract model residual standard deviation
pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())
pooled_sd

#means, SD, 95% CI per group
group_stats <- dat_all %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    `95% CI lower` = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    `95% CI upper` = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(Mean)
group_stats

#tukey HSD for pairwise comparisons
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_result$groups  # This shows letters like A, B, etc.

#using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])


SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# print all the values
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)

##Between just the mean climate variables####

dat_mean <- dat_all %>%
  filter(str_ends(ClimateVar, "_mean"))

anova_model <- aov(logR ~ ClimateVar, data = dat_mean)
summary(anova_model)

#means, SD, 95% CI per group
group_stats <- dat_mean %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    `95% CI lower` = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    `95% CI upper` = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(Mean)
group_stats

# tukey HSD for pairwise comparisons
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_groups<-tukey_result$groups  # This shows letters like A, B, etc.

# using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])
# extract ANOVA table values
anova_table <- anova_summary[[1]]
F_value <- anova_table["ClimateVar","F value"]
p_value <- anova_table["ClimateVar","Pr(>F)"]

SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# print all the values
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)

##for just growing season
dat_mean <- dat_all %>%
  filter(str_starts(ClimateVar, "GS"))

anova_model <- aov(logR ~ ClimateVar, data = dat_mean)
summary(anova_model)

# means, SD, 95% CI per group
group_stats <- dat_mean %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    `95% CI lower` = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    `95% CI upper` = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(Mean)
group_stats

# tukey HSD for pairwise comparisons
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_groups<-tukey_result$groups  # This shows letters like A, B, etc.

# using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])
# extract ANOVA table values
anova_table <- anova_summary[[1]]
F_value <- anova_table["ClimateVar","F value"]
p_value <- anova_table["ClimateVar","Pr(>F)"]

SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# print all the values
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)


###GSAI and AI_mean####
# Replace with your actual dataset
dat_sub <- dat_all %>%
  filter(ClimateVar %in% c("GSAI", "AI_mean"))

# Make factor
dat_sub$ClimateVar <- factor(dat_sub$ClimateVar, levels = c("GSAI", "AI_mean"))

# 1. One-way ANOVA
anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
anova_summary <- summary(anova_model)

# Extract ANOVA table values
anova_table <- anova_summary[[1]]
F_value <- anova_table["ClimateVar","F value"]
p_value <- anova_table["ClimateVar","Pr(>F)"]

# Pooled SD
pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())

# 2. Means, SD, 95% CI per group
group_stats <- dat_sub %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    CI_lower = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    CI_upper = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(desc(Mean))  # arrange for Tukey letters

# 3. Tukey HSD
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_groups <- tukey_result$groups %>%
  rownames_to_column("ClimateVar") 

# Using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# Or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])


SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# Print
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)


##GSPPT AND MAP####

dat_sub <- dat_all %>%
  filter(ClimateVar %in% c("GSppt", "MAP_mean"))

# Make factor
dat_sub$ClimateVar <- factor(dat_sub$ClimateVar, levels = c("GSppt", "MAP_mean"))

# 1. One-way ANOVA
anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
anova_summary <- summary(anova_model)

# Extract ANOVA table values
anova_table <- anova_summary[[1]]
F_value <- anova_table["ClimateVar","F value"]
p_value <- anova_table["ClimateVar","Pr(>F)"]

# Pooled SD
pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())

# 2. Means, SD, 95% CI per group
group_stats <- dat_sub %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    CI_lower = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    CI_upper = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(desc(Mean))  # arrange for Tukey letters

# Tukey HSD
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_groups <- tukey_result$groups %>%
  rownames_to_column("ClimateVar")

# Using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# Or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])


SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# Print
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)


##MAT and GST####


dat_sub <- dat_all %>%
  filter(ClimateVar %in% c("GStavg", "MAT_mean"))

# Make factor
dat_sub$ClimateVar <- factor(dat_sub$ClimateVar, levels = c("GStavg", "MAT_mean"))

# 1. One-way ANOVA
anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
anova_summary <- summary(anova_model)

# Extract ANOVA table values
anova_table <- anova_summary[[1]]
F_value <- anova_table["ClimateVar","F value"]
p_value <- anova_table["ClimateVar","Pr(>F)"]

# Pooled SD
pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())

# 2. Means, SD, 95% CI per group
group_stats <- dat_sub %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    CI_lower = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    CI_upper = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(desc(Mean))  # arrange for Tukey letters

# 3. Tukey HSD
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_groups <- tukey_result$groups %>%
  rownames_to_column("ClimateVar")

# Using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# Or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])


SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# Print
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)

##PET AND GSPET####

dat_sub <- dat_all %>%
  filter(ClimateVar %in% c("GSPET", "PET_mean"))

# Make factor
dat_sub$ClimateVar <- factor(dat_sub$ClimateVar, levels = c("GSPET", "PET_mean"))

# one-way ANOVA
anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
anova_summary <- summary(anova_model)

#extract ANOVA table values
anova_table <- anova_summary[[1]]
F_value <- anova_table["ClimateVar","F value"]
p_value <- anova_table["ClimateVar","Pr(>F)"]

# means, SD, 95% CI per group
group_stats <- dat_sub %>%
  group_by(ClimateVar) %>%
  summarise(
    N = n(),
    Mean = mean(logR, na.rm = TRUE),
    StDev = sd(logR, na.rm = TRUE),
    CI_lower = Mean - qt(0.975, N-1)*StDev/sqrt(N),
    CI_upper = Mean + qt(0.975, N-1)*StDev/sqrt(N)
  ) %>%
  arrange(desc(Mean))  # arrange for Tukey letters

# Tukey HSD
tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
tukey_groups <- tukey_result$groups %>%
  rownames_to_column("ClimateVar")

# using residuals
pooled_sd <- sqrt(var(residuals(anova_model)))

# or using ANOVA table (more standard)
anova_table <- summary(anova_model)[[1]]
pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])


SS_total <- sum(anova_table[,"Sum Sq"])
SS_between <- anova_table["ClimateVar", "Sum Sq"]
R_sq <- SS_between / SS_total
R <- sqrt(R_sq)

# Print all the values
list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
     R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)


##---END OF CODE---####