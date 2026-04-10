##Created by Nidhi Vinod
##Priorities of woody species trait-climate associations at continental scale
##Updated on April,9, 2026
##Contains the entire analysis for the MS for the site climate data
#------------------------------------------------------------------------------------------------- 
##load libraries
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

  #------------------------------------------------------------------------------------------------- 
  
  #set directory      
  #------------------------------------------------------------------------------------------------- 
  
  #load data
  
  site<-read.csv("~/Desktop/Traits-and-Climate/Raw data/archive/NEON-CTFS_02May25_sites.csv")
  site_climate<-read.csv("~/Desktop/Traits-and-Climate/Raw data/species_16Feb26.csv")
  
  #delete empty row
  site_climate<-site_climate[-405,]
  
  #such a huge dataset that I am putting the colnames into a df to see which to remove
  # names<-colnames(site_climate)
  # names<-as.data.frame(names)
  
  site<-site[,c(1,2,12,14,19,20,21,24,41,42)]
  
  #filtering data as suggested by Lawren
  new_dat<-site_climate[,c(1:4,35,41,44,47,65:71,73:76,78:80,
                           82,83,89,94:117,119:121,124,125,131,
                           134,136,138,147:149,156,157)]
  
  site<-site%>%rename("AI_mean"="ai", "PET_mean"="pet", "MAT_mean" = "MAT", "MAP_mean"="MAP")
  
  #write.csv(new_dat,"~/Desktop/Traits-and-Climate/Raw data/trait_climate_gs_dat.csv")
  
  ##replace the species MAT and MAP and other climatw with data of the site climate
  
  # Define the climate variable names
  climate_vars <- c("GSAI", "GSPET", "GSppt", "GStavg", "AI_mean", "PET_mean", "MAT_mean", "MAP_mean")
  
  # Replace the 8 climate variables in `site_climate` with values from `site` based on site
  new_dat[climate_vars] <- site[match(new_dat$Site_code, site$Site_code), climate_vars]
  
  #remove white mountain sites because it mostly has herbacious species
  dat_clean <- new_dat %>% 
    filter(!Site %in% c("WhiteMountains", "Laupahoehoe", "Palamanui"))
  
  
  #write.csv(dat_clean,"species_clim_changed_to_site_clim.csv")
 
   #filter out site and site code for cormat####
  data<-dat_clean[,-c(1:4)]
  
  ##Correlation matrix:
  
  setwd("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii")
  
  data<-as.data.frame(data)#convert to dataframe
  data<-abs(data)
  
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
  # corData_1<-lapply(corData,as.numeric)
  # corData_1 <- data.frame(matrix(unlist(corData_1), nrow=length(corData_1), byrow=FALSE))
  # colnames(corData_1) <- colnames(corData)
  # corData<-corData_1
  # pData_1<-lapply(pData,as.numeric)
  # pData_1 <- data.frame(matrix(unlist(pData_1), nrow=length(pData_1), byrow=FALSE))
  # colnames(pData_1) <- colnames(pData)
  # pData<-pData_1
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
  
  # corData_1<-lapply(corData,as.numeric)
  # corData_1 <- data.frame(matrix(unlist(corData_1), nrow=length(corData_1), byrow=FALSE))
  # colnames(corData_1) <- colnames(corData)
  # corData<-corData_1
  # pData_1<-lapply(pData,as.numeric)
  # pData_1 <- data.frame(matrix(unlist(pData_1), nrow=length(pData_1), byrow=FALSE))
  # colnames(pData_1) <- colnames(pData)
  # pData<-pData_1
  
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
  
  ###Get Highest R values for significant correaltions for Table S_10####
  
  #load data
  raw<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RAWcor.csv")
  rank<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RANKcor.csv")
  log<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/LOGcor.csv")
  
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
  
  #final_mat<-pmax(mat_1, mat_2, mat_3)
  
  max_matrix <- matrix(0, nrow = 67, ncol = 68)
  
  # Loop through each element to find the maximum
  for (i in 1:nrow(max_matrix)) {
    for (j in 1:ncol(max_matrix)) {
      max_matrix[i, j] <- max(a[i, j], b[i, j], c[i, j])
    }
  }
  
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
  
  
  #write.csv(max_matrix,"~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/Rmax_values_si_feb26.csv")
  
  ##for calculating the length of traits and such
  
  raw<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RAWcor_sig.csv")
  rank<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RANKcor_sig.csv")
  log<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/LOGcor_sig.csv")
  rawp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RAWp.csv")
  rankp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RANKp.csv")
  logp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/LOGp.csv")

  
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
  
  
  #it doesn't matter which df we use to add the high r value
  raw_max$High_R <- pmax( #gets the highest R values
    raw_lon$Value,
    log_lon$Value,
    rank_lon$Value,
    na.rm = TRUE
  )
  colnames(raw_max)[4] <- "Srl"
  
  #calculate the number of sig correlations
  corr_long<-raw_max[,-c(3,4)]
  
  climate_vars <- c("GSAI", "GSPET",
                    "GSppt", "GStavg",
                    "AI_mean", "PET_mean",
                    "MAT_mean", "MAP_mean")
  
  # Helper: not-in operator
  `%notin%` <- Negate(`%in%`)
  
  # Keep only trait × climate correlations
  
  corr_tc <- corr_long %>%
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
  
  # 2) Keep only plant traits
  plant_traits <- unique(corr_long$Row[!(corr_long$Row %in% climate_vars)])
  
  # 3) Count plant traits
  length(plant_traits)  #
  
  # 6) Numbers for manuscript
  
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
  
  # Define lists of trait and climate variable names
  trait_vars <- trait_vars <- c(
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
    "SeedMassTRY_mg")
  
  climate_vars <- c("GSAI", "GSPET", "GSppt", "GStavg", "AI_mean", "PET_mean", "MAT_mean", "MAP_mean")
  
  # raw_max_clean <- raw_max %>%
  #   # Keep only climate × trait correlations
  #   filter(
  #     Col %in% climate_vars & Row %in% trait_vars
  #   )
  raw_max_clean<-raw_max[-c(1:3953),]
  
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
  
  write.csv(raw_max_n,"~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/highest_r_for_site_sp_metawin.csv")
  
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
  
  #calculating yes and no for every cliamte variable
  yes_count <- data_sorted_sp %>%
    group_by(S_A) %>%
    summarise(Yes_Count = sum(y_n == "YES"))
  
  no_count <- data_sorted_sp %>%
    group_by(S_A) %>%
    summarise(No_Count = sum(y_n == "NO"))
  
  yes_no<-cbind(yes_count,no_count$No_Count)
  
  yes_no<-yes_no %>% rename("No_Count" = "no_count$No_Count")
  
  yes_no$total<-yes_no$Yes_Count+yes_no$No_Count
  
  yes_no$perc_sig<-yes_no$Yes_Count/yes_no$total
  
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
  
  ##Calculate proportion of Significance#####
  
  dat<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/HighlightTable.csv")
  
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
  
  # #run the function for every data set
  dat_lon<-df(dat)# converting matrix to rows, cols from fn above
  
  
  #name the fourth col
  colnames(dat_lon)[4] <- "Srl"
  
  #removing trait-trait correlations
  #follow the second column till you see gsppt, cut until there
  # Define lists of trait and climate variable names
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
    "SeedMassTRY_mg")
  
  
  
  climate_vars <- c("GSAI", "GSPET", "GSppt", "GStavg", "AI_mean", "PET_mean", "MAT_mean", "MAP_mean")
  
  dat_lon<-dat_lon[-c(1:3953),]
  
  #remove climate-climate correlations and remove those numbers in every 
  #repetition
  dat_lon_clean<- dat_lon[!(raw_max$Srl %in% 60:67), ]
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
  
  #calculating yes and no for every cliamte variable
  yes_count <- dat_lon_clean %>%
    group_by(S_A) %>%
    summarise(Yes_Count = sum(Value == "YES"))
  
  no_count <- dat_lon_clean %>%
    group_by(S_A) %>%
    summarise(No_Count = sum(Value == "NO"))
  
  yes_no<-cbind(yes_count,no_count$No_Count)
  
  yes_no<-yes_no %>% rename("No_Count" = "no_count$No_Count")
  
  yes_no$total<-yes_no$Yes_Count+yes_no$No_Count
  
  yes_no$perc_sig<-yes_no$Yes_Count/yes_no$total
  
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
  write.csv(results,"~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/results_for_swtiching_siteclimate_for_sp_wo_hawaii.csv")
  #this above analysis is the right one.
  
  #Metawin Analysis#####
  
  #r2_table<-read.csv("~/Desktop/Traits-and-Climate/Clean data/Nov_5_species_climate_to_site/highest_r_for_site_sp_metawin.csv")#adding the r2 table from meta-win into a new df
  
  raw<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RAWcor_sig.csv")
  rank<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RANKcor_sig.csv")
  log<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/LOGcor_sig.csv")
  rawp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RAWp.csv")
  rankp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/RANKp.csv")
  logp<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/LOGp.csv")
  
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
  #write.csv(r2_table, "~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/r_table_clean_site_clim.csv")
  
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
  
  #r2_table$R_abs <- abs(r2_table$R_Value)
  
  
  pairwise <- t(combn(climate_vars, 2))
  logrr_data <- data.frame()
  
  for (i in 1:nrow(pairwise)) {
    var1 <- pairwise[i, 1]
    var2 <- pairwise[i, 2]
    for (trait in trait_vars) {
      r2_1 <- r2_table$R[r2_table$Trait == trait & r2_table$Climate_Variable == var1]
      r2_2 <- r2_table$R[r2_table$Trait == trait & r2_table$Climate_Variable == var2]
      
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
  
  
  #Meta-analysis using metafor
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
  write.csv(my_r_data,"~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/meta_win_summary_high_r_data.csv")
  
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
  write.csv(summary_table_p_my_data,"~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/metawin_siteclimate_for_sp_wo_hawaii.csv")

    # Forest plot for first comparison
  forest(results[[1]], main = names(results)[1])
  ##############################
  
  #separating the R values so I can use it for the Table S_13 in the paper

  # Calculate mean R for each climate variable pair
  dat<-r2_table
  mean_Rs <- dat %>%
    group_by(Climate_Variable) %>%   # your climate variable column names here
    summarise(mean_R = mean(abs(R_Value), na.rm = TRUE)) %>%
    ungroup()
  
  print(mean_Rs)
  write.csv(mean_Rs,"mean_r_values.csv")
  
  ##Making the final table:#####
  
  library(dplyr)
  library(tidyr)
  
  # Load your previously created objects ---
  # results -> prop test results
  # summary_table_p_my_data -> LogRR + CI + p-values
  # mean_Rs -> mean R for each climate variable
  
  # First, separate Climate_Variable_1 and Climate_Variable_2 from Comparison
  logrr_df <- summary_table_p_my_data %>%
    separate(Comparison, into = c("Climate_Variable_1", "Climate_Variable_2"), sep = "_vs_")
  
  # Merge prop-test results with LogRR results
  combined <- results %>%
    separate(Comparison, into = c("Climate_Variable_1", "Climate_Variable_2"), sep = " vs ") %>%
    left_join(logrr_df, by = c("Climate_Variable_1", "Climate_Variable_2"))
  
  # Mapping between names in pairwise/comparison and mean_Rs
  name_map <- data.frame(
    Combined = c("GSP", "GST", "GSPET", "GSAI", "MAAI", "MAPET", "MAT", "MAP"),
    MeanR = c("GSppt", "GStavg", "GSPET", "GSAI", "AI_mean", "PET_mean", "MAT_mean", "MAP_mean")
  )
  
  # Merge Mean_R1
  combined <- combined %>%
    left_join(name_map, by = c("Climate_Variable_1" = "Combined")) %>%
    left_join(mean_Rs, by = c("MeanR" = "Climate_Variable")) %>%
    rename(Mean_R1 = mean_R) %>%
    dplyr::select(-MeanR)
  
  # Merge Mean_R2
  combined <- combined %>%
    left_join(name_map, by = c("Climate_Variable_2" = "Combined")) %>%
    left_join(mean_Rs, by = c("MeanR" = "Climate_Variable")) %>%
    rename(Mean_R2 = mean_R) %>%
    dplyr::select(-MeanR)
  
  # Round numeric columns for nicer display
  combined <- combined %>%
    mutate(across(c(Proportion1, Proportion2, p_value, LogRR_Estimate, CI_Lower, CI_Upper, P_Value, Mean_R1, Mean_R2),
                  ~round(., 3)))
  
  #  Reorder columns
  combined <- combined %>%
    dplyr::select(Climate_Variable_1, Climate_Variable_2, Proportion1, Proportion2, p_value, Significance,
                  Mean_R1, Mean_R2, LogRR_Estimate, CI_Lower, CI_Upper, P_Value, Interpretation)
  
  # View final table
  print(combined)
  
  
  # Separate prop-test results
  results_sep <- results %>%
    separate(Comparison, into = c("Climate_Variable_1", "Climate_Variable_2"), sep = " vs ")
  
  logrr_df <- summary_table_p_my_data %>%
    separate(Comparison, into = c("Climate_Variable_1", "Climate_Variable_2"), sep = "_vs_")
  
  
  # Convert long names in logrr_df to short names
  logrr_df_short <- logrr_df %>%
    mutate(
      Climate_Variable_1 = case_when(
        Climate_Variable_1 == "GSppt" ~ "GSP",
        Climate_Variable_1 == "GStavg" ~ "GST",
        Climate_Variable_1 == "GSPET" ~ "GSPET",
        Climate_Variable_1 == "GSAI" ~ "GSAI",
        Climate_Variable_1 == "AI_mean" ~ "MAAI",
        Climate_Variable_1 == "PET_mean" ~ "MAPET",
        Climate_Variable_1 == "MAT_mean" ~ "MAT",
        Climate_Variable_1 == "MAP_mean" ~ "MAP",
        TRUE ~ Climate_Variable_1
      ),
      Climate_Variable_2 = case_when(
        Climate_Variable_2 == "GSppt" ~ "GSP",
        Climate_Variable_2 == "GStavg" ~ "GST",
        Climate_Variable_2 == "GSPET" ~ "GSPET",
        Climate_Variable_2 == "GSAI" ~ "GSAI",
        Climate_Variable_2 == "AI_mean" ~ "MAAI",
        Climate_Variable_2 == "PET_mean" ~ "MAPET",
        Climate_Variable_2 == "MAT_mean" ~ "MAT",
        Climate_Variable_2 == "MAP_mean" ~ "MAP",
        TRUE ~ Climate_Variable_2
      )
    )
  
  # Now join with results using short n
  
  # Join with LogRR + significance table
  combined <- results_sep %>%
    left_join(logrr_df_short, by = c("Climate_Variable_1", "Climate_Variable_2"))
  
  # Merge Mean_R1 and Mean_R2 using case_when
  combined <- combined %>%
    mutate(
      Mean_R1 = case_when(
        Climate_Variable_1 == "GSP" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GSppt"],
        Climate_Variable_1 == "GST" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GStavg"],
        Climate_Variable_1 == "GSPET" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GSPET"],
        Climate_Variable_1 == "GSAI" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GSAI"],
        Climate_Variable_1 == "MAAI" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "AI_mean"],
        Climate_Variable_1 == "MAPET" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "PET_mean"],
        Climate_Variable_1 == "MAT" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "MAT_mean"],
        Climate_Variable_1 == "MAP" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "MAP_mean"],
        TRUE ~ NA_real_
      ),
      Mean_R2 = case_when(
        Climate_Variable_2 == "GSP" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GSppt"],
        Climate_Variable_2 == "GST" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GStavg"],
        Climate_Variable_2 == "GSPET" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GSPET"],
        Climate_Variable_2 == "GSAI" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "GSAI"],
        Climate_Variable_2 == "MAAI" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "AI_mean"],
        Climate_Variable_2 == "MAPET" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "PET_mean"],
        Climate_Variable_2 == "MAT" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "MAT_mean"],
        Climate_Variable_2 == "MAP" ~ mean_Rs$mean_R[mean_Rs$Climate_Variable == "MAP_mean"],
        TRUE ~ NA_real_
      )
    )
  
  
  
  # Round numeric columns
  combined <- combined %>%
    mutate(across(c(Proportion1, Proportion2, p_value, LogRR_Estimate, CI_Lower, CI_Upper, P_Value, Mean_R1, Mean_R2),
                  ~round(., 3)))
  
  # Reorder columns
  combined <- combined %>%
    dplyr::select(Climate_Variable_1, Climate_Variable_2, Proportion1, Proportion2, p_value, Significance,
                  Mean_R1, Mean_R2, LogRR_Estimate, CI_Lower, CI_Upper, P_Value, Interpretation)
  
  combined <- combined %>%
    mutate(
      LogRR_Estimate = case_when(
        Climate_Variable_1 == "MAP" & Climate_Variable_2 == "MAT" ~ -0.258,
        TRUE ~ LogRR_Estimate
      ),
      CI_Lower = case_when(
        Climate_Variable_1 == "MAP" & Climate_Variable_2 == "MAT" ~ -0.565,
        TRUE ~ CI_Lower
      ),
      CI_Upper = case_when(
        Climate_Variable_1 == "MAP" & Climate_Variable_2 == "MAT" ~ 0.049,
        TRUE ~ CI_Upper
      ),
      P_Value = case_when(
        Climate_Variable_1 == "MAP" & Climate_Variable_2 == "MAT" ~ 0.0995,
        TRUE ~ P_Value
      ),
      Interpretation = case_when(
        Climate_Variable_1 == "MAP" & Climate_Variable_2 == "MAT" ~ "No significant difference (p ≥ 0.05)",
        TRUE ~ Interpretation
      )
    )
  
  #reorder like in my word doc:
  library(dplyr)
  
  # Your desired order as a list of pairs
  desired_order <- list(
    c("GSP", "GSPET"),
    c("GSP", "GST"),
    c("GSP", "GSAI"),
    c("GSPET", "GSAI"),
    c("GSPET", "GST"),
    c("GST", "GSAI"),
    
    c("MAAI", "MAP"),
    c("MAPET", "MAAI"),
    c("MAAI", "MAT"),
    c("MAP", "MAT"),
    c("MAPET", "MAP"),
    c("MAPET", "MAT"),
    
    c("MAAI", "GSAI"),
    c("MAP", "GSP"),
    c("MAPET", "GSPET"),
    c("GST", "MAT")
  )
  
  # Create a helper column that is always sorted pair
  combined <- combined %>%
    rowwise() %>%
    mutate(
      Comparison_pair = list(c(Climate_Variable_1, Climate_Variable_2))
    ) %>%
    mutate(
      Comparison_index = match(TRUE, sapply(desired_order, function(x) all(x == Comparison_pair)))
    ) %>%
    ungroup() %>%
    arrange(Comparison_index) %>%
    dplyr::select(-Comparison_pair, -Comparison_index)
  
  write.csv(combined,"~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/table3_for_sp_clim_replaced_with_site_clim_wo_hawaii.csv")
  
##Mean,SE, SD####
  
  dat<-read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/highest_r_for_site_sp_metawin.csv")
  
  #remove rows that have NO
  yes_s<- data_sorted_sp %>% 
    filter(!y_n %in% "NO")
  
  
  mean_sd_se <- yes_s %>%
    group_by(S_A) %>%  # your grouping variable
    summarise(
      mean_r = mean(abs(Number), na.rm = TRUE),
      sd_r = sd(abs(Number), na.rm = TRUE),
      se_r = sd(abs(Number), na.rm = TRUE) / sqrt(sum(!is.na(Number)))
    ) %>%
    ungroup()
  
  write.csv(mean_sd_se, "~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/mean_se_sd_for_site_sp_climate_replaced.csv")
  
  #### Multiple regressions: predicting traits from MAP and MAT ####
  #### Camila D. Medeiros ####
  #### 21 Oct 25 ####
  
  library(car)
  library(MASS)
  library(MuMIn) #package used to calculate AICc
  library(hier.part)
  library(rdacca.hp)
  
  # Set WD:
  
  # Load data:
  data<-dat_clean
  #data <- read.csv("~/Desktop/Traits-and-Climate/Raw data/trait_climate_gs_dat.csv")
  #data<-data[,-1] # trait data
  
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
  
  
  
  #### TLP ####
  ##### a) Fit linear model usin OLS: #####
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
  
  ##### b) Hierarchical partitioning (to find out the percentage of contribution from @ predictor variable): #####
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

  ###Nidhi's modified code####
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
  write.csv(RAW.df, "MultipleRegressionRes_raw_for_species_site_wo_hawaii.csv", row.names = FALSE)
  
  ####Nidhi's modified code####
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
  write.csv(LOG.df, "MultipleRegressionRes_log_sp_site.csv", row.names = FALSE)
  
  #setwd("~/Desktop/Traits-and-Climate/Clean data/without hawaii sites analysis")
  ##now find the best fit:
  ### 1. Load the two CSVs
  raw <- read.csv("MultipleRegressionRes_raw_for_species_site_wo_hawaii.csv")
  log <- read.csv("MultipleRegressionRes_log_sp_site.csv")
  
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
  
  ### 7. Write final table
  write.csv(merged, "Merged_Raw_Log_Comparison.csv", row.names = FALSE)
  
  # ### 8. Print summary
  # cat("Summary:\n")
  # cat("Mean % contribution MAP:", round(mean_MAP_contrib, 2), "\n")
  # cat("Mean % contribution MAT:", round(mean_MAT_contrib, 2), "\n")
  # cat("% traits with MAP contribution higher:", round(perc_MAP_higher*100, 2), "%\n")
  # cat("% traits with MAT contribution higher:", round(perc_MAT_higher*100, 2), "%\n")
  # 
  ### 9. Optional: save summary stats as a small table
  summary_df <- data.frame(
    Summary = c("mean percent contribution MAP", "mean percent contribution MAT",
                "% traits with MAP contribution higher", "% traits with MAT contribution higher",
                "SE % contribution MAT:", "SE % contribution MAP:"),
    Value = c(mean_MAP_contrib, mean_MAT_contrib, perc_MAP_higher*100, perc_MAT_higher*100,se_MAT,
              se_MAP )
  )
  write.csv(summary_df, "Merged_Raw_Log_Summary.csv", row.names = FALSE)
  
  ## Highest R cvalue for every trait df####
  
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
  
  setwd("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii")
  
  write.csv(correlation_counts,"sp_site_without_hawaii_n_correlations.csv")
  write.csv(hr_sep,"sp_site_without_hawaii_r_values.csv")
  
  #Figure 3#### R mean values with SE####
  
  library(dplyr)
  library(ggplot2)
  
  r_dat <- read.csv("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii/mean_se_sd_for_site_sp_climate_replaced.csv")
  
  r_dat<-r_dat%>%rename ("Climate" = "S_A")
  # Add MA vs GS grouping and CLD letters
  r_dat <- r_dat %>%
    mutate(
      GS_MA = ifelse(grepl("^GS", Climate), "GS", "MA"),
      CLD = case_when(
        Climate == "MAAI"  ~ "A",
        Climate == "MAP"   ~ "A",
        Climate == "MAPET" ~ "A",
        Climate == "MAT"   ~ "B",
        Climate == "GSP"   ~ "ab",
        Climate == "GST"   ~ "ab",
        Climate == "GSPET" ~ "a",
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
  
  # Set1 + orange palette 
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
  
  #pie chart
  
  values <- c(MAP = 79.6610169491525, MAT = 20.3389830508475)
  
  colors <- c(
    "MAP" = "#4575B4", # blue
    "MAT" = "#D73027")  # red)  
  
  pie3D(values,
        labels = NA,          # no names or values
        explode =  0.1, # slightly separate slices
        theta = 1.5,          # faces more toward screen
        col = colors)
  #------------------------------------------------------------------------------------------------- 
  ##CV, range, slrt and leven's test for MAT and MAP####
  
  #load data
  
  # dat_clean generated above
  
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
  
  
  ##ANOVA for Table S_14#####
  library(tidyverse)
  library(car)
  library(emmeans)
  library(multcompView)
  library(agricolae)
  
 #set working directory
  setwd("~/Desktop/Traits-and-Climate/Clean data/sp_site_without_hawaii")
  
  #load data
  df <- read.csv("sp_site_without_hawaii_r_values.csv")
  
  # Log-transform r-values if needed

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
  
  # log-transform r-values
  dat_all <- dat %>% mutate(logR = log(R) + 1e-6)
  
  # One-way ANOVA
  anova_model <- aov(logR ~ ClimateVar, data = dat_all)
  
  #ANOVA summary
  anova_summary <- summary(anova_model)
  anova_summary
  
  # Extract model residual standard deviation
  pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())
  pooled_sd
  
  # Means, SD, 95% CI per group
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
  
  # 5. Tukey HSD for pairwise comparisons
  tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
  tukey_result$groups  # This shows letters like A, B, etc.
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
  
##Between just Mean climate variables####
  
  dat_mean <- dat_all %>%
    filter(str_ends(ClimateVar, "_mean"))
  
  anova_model <- aov(logR ~ ClimateVar, data = dat_mean)
  summary(anova_model)
  
  # Means, SD, 95% CI per group
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
  
  # Tukey HSD for pairwise comparisons
  tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
  tukey_groups<-tukey_result$groups  # This shows letters like A, B, etc.
  
  # Using residuals
  pooled_sd <- sqrt(var(residuals(anova_model)))
  
  # Or using ANOVA table (more standard)
  anova_table <- summary(anova_model)[[1]]
  pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])
  # Extract ANOVA table values
  anova_table <- anova_summary[[1]]
  F_value <- anova_table["ClimateVar","F value"]
  p_value <- anova_table["ClimateVar","Pr(>F)"]
  
  SS_total <- sum(anova_table[,"Sum Sq"])
  SS_between <- anova_table["ClimateVar", "Sum Sq"]
  R_sq <- SS_between / SS_total
  R <- sqrt(R_sq)
  
  # Print
  list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
       R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)
  
  ##For just growing season variables####
  dat_mean <- dat_all %>%
    filter(str_starts(ClimateVar, "GS"))
  
  anova_model <- aov(logR ~ ClimateVar, data = dat_mean)
  summary(anova_model)
  
  #Means, SD, 95% CI per group
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
  
  #Tukey HSD for pairwise comparisons
  tukey_result <- HSD.test(anova_model, "ClimateVar", group = TRUE)
  tukey_groups<-tukey_result$groups  # This shows letters like A, B, etc.
  
  # Using residuals
  pooled_sd <- sqrt(var(residuals(anova_model)))
  
  # Or using ANOVA table (more standard)
  anova_table <- summary(anova_model)[[1]]
  pooled_sd <- sqrt(anova_table["Residuals", "Mean Sq"])
  # Extract ANOVA table values
  anova_table <- anova_summary[[1]]
  F_value <- anova_table["ClimateVar","F value"]
  p_value <- anova_table["ClimateVar","Pr(>F)"]
  
  SS_total <- sum(anova_table[,"Sum Sq"])
  SS_between <- anova_table["ClimateVar", "Sum Sq"]
  R_sq <- SS_between / SS_total
  R <- sqrt(R_sq)
  
  # Print
  list(F_value = F_value, p_value = p_value, pooled_sd = pooled_sd,
       R_sq = R_sq, R = R, group_stats = group_stats, tukey_groups = tukey_groups)
  
  
  ###Numeric logR correlations for GSAI and AI_mean####
  # Replace with your actual dataset
  dat_sub <- dat_all %>%
    filter(ClimateVar %in% c("GSAI", "AI_mean"))
  
  # Make factor
  dat_sub$ClimateVar <- factor(dat_sub$ClimateVar, levels = c("GSAI", "AI_mean"))
  
  #One-way ANOVA
  anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
  anova_summary <- summary(anova_model)
  
  # Extract ANOVA table values
  anova_table <- anova_summary[[1]]
  F_value <- anova_table["ClimateVar","F value"]
  p_value <- anova_table["ClimateVar","Pr(>F)"]
  
  # Pooled SD
  pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())
  
  #Means, SD, 95% CI per group
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
  
  #Tukey HSD
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
  
  
  ##MAP and GSP####
  
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
  
  # Means, SD, 95% CI per group
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
  
  #Tukey HSD
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
  
  # One-way ANOVA
  anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
  anova_summary <- summary(anova_model)
  
  # Extract ANOVA table values
  anova_table <- anova_summary[[1]]
  F_value <- anova_table["ClimateVar","F value"]
  p_value <- anova_table["ClimateVar","Pr(>F)"]
  
  # Pooled SD
  pooled_sd <- sqrt(anova_model[[1]][["Residuals"]] %>% var())
  
  # Means, SD, 95% CI per group
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
  
  ##PET AND GSPET####
  
  dat_sub <- dat_all %>%
    filter(ClimateVar %in% c("GSPET", "PET_mean"))
  
  # Make factor
  dat_sub$ClimateVar <- factor(dat_sub$ClimateVar, levels = c("GSPET", "PET_mean"))
  
  #One-way ANOVA
  anova_model <- aov(logR ~ ClimateVar, data = dat_sub)
  anova_summary <- summary(anova_model)
  
  # Extract ANOVA table values
  anova_table <- anova_summary[[1]]
  F_value <- anova_table["ClimateVar","F value"]
  p_value <- anova_table["ClimateVar","Pr(>F)"]
  
  #Means, SD, 95% CI per group
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
  
  #Tukey HSD
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
  
  