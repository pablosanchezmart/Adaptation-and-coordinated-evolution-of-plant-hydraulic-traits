################################################## MERGE METHODOLOGY ###############################################################

# Sanchez Pablo

### PACKAGES ####
library(sp)
library(stringr)
library(ALA4R)
library(rgbif)
library(raster)
library(futile.logger)
library(doBy)
#remove.packages("speciesmap")
#devtools::install_github("remkoduursma/speciesmap")
library(speciesmap)
#https://remkoduursma.github.io/speciesmap/articles/Using_speciesmap.html
library(knitr)
library(dplyr)
library(dismo)
library(plyr)

setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")

### FUNCTIONS ####

## Mode function (for numerical variables aggregation)
mean.fun <- function(x){
  y <- mean(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}
## Mode function (for categorical variables aggregation)
mode.fun <- function(v) {
  v <- as.character(v)
  uniqv <- unique(v)
  mode <- uniqv[which.max(tabulate(match(v, uniqv)))]
  if(is.na(mode)){                                            # To avoid NAs
    uniqv <- uniqv[which(!is.na(uniqv))]
    mode <- uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  return(mode)
}

# To combine already existing variables
combFun <- function(df, vrbls){
  for(variable in vrbls){
    cat("Scrapping", variable, "\n")
    varX <- paste(variable, ".x", sep = "")
    varY <- paste(variable, ".y", sep = "")
    for(i in 1:length(df$Species))
      if(is.na(df[i , varX])){
        df[i , variable] <- df[i , varY]
      } else {
        df[i , variable] <- df[i , varX]
      }
  } # end creating new columns
  return(df)
} # end function

################################ MERGE HYDRATRY AND LIU ET AL. 2019 ################################################################

### PREPRING HYDRATRY AND LIU ET AL. 2019 DATASETS TO LET THEM MERGE ####

## Charge HydraTRY database...
ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T) # Without liu species
names(ht)
# Select the non-redundant variables (without the logaritmic transformations)
ht <- ht[, -c(1, 60:74)]

## Charge Liu et al. 2019 database...
liuData <- readxl::read_excel("Others/Liu_etal-2019_Data.xlsx", sheet = "full_data")
#liuData <- read.csv("Liu_etal-2019_Data_original.csv", header = T)                
names(liuData)
# Select variables
liuData <- liuData[, -c(1, 2, 9, 10, 11, 31, 33)]
liuData$Species <- as.factor(liuData$Species)

### Preparation of database (same names)

names(liuData)[names(liuData) == "Ppd"] <- "MinWP_pd"
names(liuData)[names(liuData) == "Pmd"] <- "MinWP_md"
names(liuData)[names(liuData) == "TLP"] <- "PItlp"
names(liuData)[names(liuData) == "Life form"] <- "Life.form"
names(liuData)[names(liuData) == "Leaf form"] <- "Leaf.form"
names(liuData)[names(liuData) == "Biome"] <- "biome"
names(liuData)[names(liuData) == "Group"] <- "group"
names(liuData)[names(liuData) == "Family"] <- "family"
names(liuData)[names(liuData) == "Genus"] <- "genus"

# Note that Kl needs to be divided by 10^4 to have same units (Liu_etal-2019_Data.csv already changed)
liuData$Kl <- liuData$Kl/(10^4)

# Change level names for biome factor
liuData$biome <- plyr::revalue(liuData$biome, c("BOR"= "Boreal forest", "DES" = "Desert", "TMR" = "Temperate rain forest", "TMS" = "Temperate forest", "TRR" = "Tropical rain forest",
                                                "TRS" = "Tropical seasonal rainforest","WDS" = "Woodland/shrubland"))
# Change level names for group factor
liuData$group <- plyr::revalue(liuData$group, c("Ang" = "Angiosperms", "Gym" = "Gymnosperms")) # Only this two levels

# Hv calculation for Liu
liuData$Hv <- 1/liuData$AlAs

### Aggregate Liu dataset
#Numerical variables
numVar <- c("Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "PItlp", "AlAs", "WD", "Vdia", "VD", "VLmax", "ALT", "AI", "Hv")
# Categorical variables
catVar <- c("Life.form", "Leaf.form", "Soiltype", "biome", "family", "genus", "group")

liuDataNum <- aggregate(liuData[, numVar], by = list(liuData$Species), FUN = mean.fun)  # Had problems with the previus mean fun (not account for NAs omission)
liuDataCat <- aggregate(liuData[, catVar], by = list(liuData$Species), FUN = mode.fun)

# Merge categorical and numerical variables
liuData <- merge(liuDataNum, liuDataCat, by = "Group.1")

names(liuData)[names(liuData) == "Group.1"] <- "Species"
names(liuData)

### MERGE HYDRATRY AND LIU ETAL 2019####

ht_liu <- join(ht, liuData, by = "Species", type = "full")
head(ht_liu)
ht_liu[which(ht_liu$Species == "Abies alba"), ]
ht[which(ht$Species == "Abies alba"), ]

names(ht)
length(ht_liu_merge$Species)
names(ht_liu)

# Bucle to add new values to Ht dataset (x)

variables <- names(liuData[, -1])

# Create new columns with completed obs for each variable
ht_liu_merge <- combFun(ht_liu_merge, variables)

# Checking
# with(ht_liu_merge, length(which(!is.na(Hmax))))
# with(ht_liu_merge, length(which(!is.na(Kl.y))))
# length(ht_liu_merge$Species)
# new <- with(ht_liu_merge, length(which(is.na(Ks.x) & !is.na(Ks.y)))) # New values
# old <- with(ht_liu_merge, length(which(!is.na(Ks.x)))) # Old values
# new + old # Total

### SELECT COMPLETED VARIABLES AND SAVE COMBINED DATASET ####

# Select complete variables
ht_liu <- ht_liu_merge[,c("Species", "genus", "family", "order", "group", "Growth.form", "Life.form", "Leaf.form", "biome", "Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "AlAs", "WD", "Vdia", "VD", "VLmax", "ALT", "AI",
                    "Hv", "HSM", "PItlp", "WUE", "SLA", "Aarea", "Amass", "LL", "pLL", "N", "Narea", "Rd", "lat_mean", "long_mean", "MAT_mean", "MAP_mean", "tmin", "tmax", "srad", "wind", "vapr", "TMean.year", "Diurnal.T.Range", 
                    "Isothermality", "Seasonality", "TMax.warmest", "TMin.Coldest" , "TRange", "TMean.Wettest", "TMean.Driest", "TMean.Warmest", "TMean.Coldest", "Prec.year", "Prec.Wettest", "Prec.Driest", "Prec.Seasonality", "Prec.WetQ",
                    "Prec.DryQ", "Prec.WarmQ", "Prec.ColdQ", "PPET", "AbsoluteDepthBedRock_1km_ll", "PredictedClassTAXOUSDA_1km_ll",  "PredictedClassTAXNWRB_1km_ll", "Soiltype")]
length(ht_liu)
# 68
names(ht_liu)
length(ht_liu$Species)
#2225
summary(ht_liu)
write.csv(ht_liu, "HydraTRY_Pablo_logData_HmaxCompl_Liu.csv", row.names = F)


################################ ADD gs VALUES TO HT ################################################################################

ht <- read.csv("HT_Pablo_2019_log.csv", header = T)
gs <- readxl::read_excel("Others/gsDataset.xlsx")

gs$gmax <- as.numeric(gs$..16)
gs$Species <- as.factor(gs$..11)
gs$Species <- str_replace_all(gs$Species, "_", " ")
gs$Species <- word(gs$Species, 1, 2)
unique(gs$Species)

## Agregate data

gsData <- aggregate(gs[, "gmax"], by = list(gs$Species), FUN = mean.fun)
names(gsData)[names(gsData) == "Group.1"] <- "Species"

### MERGE ####
ht_gs_merged <- merge(ht, gsData, by = "Species", all = TRUE)
length(ht_gs_merged$Species) - length(ht$Species) # New species (gs)

# Checking
names(ht_gs_merged)
with(ht_gs_merged, length(which(is.na(gmax))))
length(ht_gs_merged$Species)
length(ht_gs_merged)
summary(ht_gs_merged)
summary(ht)
# Save complete dataframe
write.csv(ht_gs_merged, "HT_Pablo_2019_log_gmax.csv", row.names = F)


################################ MERGE HYDRATRY WITH PREVIOUSLY EXTRACTED CLIMATE DATA ################################################################################

ht <- read.csv("HT_Pablo_2019_log_gmax.csv", header = T)                # Last version
length(ht$Species)
clim <- read.csv("Climate_data/climate_data_extracted.csv", header = T)
summary(clim)
# Delete previous climate data
names(ht)
ht <- ht[, -c(37:84)]
names(ht)

### MERGE ###
ht_clim_merge <- merge(ht, clim, by = "Species", all = TRUE)

# Checking
names(ht_clim_merge)
with(ht_clim_merge, length(which(is.na(Hmax))))
length(ht_clim_merge$Species)
summary(ht_clim_merge$TMax.warmest)


# Save complete dataframe
write.csv(ht_clim_merge, "HT_Pablo_2019_log_gmax_clim.csv", row.names = F)
