#### MEAN LONGITUDE AND LATITUDE #############################################################################################################

# P Sanchez-Martinez

print("Calculating mean longitude and latitude for some plots...")

### PACKAGE ####

library(sp)
library(stringr)
library(ALA4R)
library(rgbif)
library(raster)
library(plyr)

### FUNCTIONS ####

meanLongLatFun <- function(species = ht$Species, occ.pth = "data/species_data/occurrences/processed_data_outl/"){
    print(paste0("extracting longitude and latitude from ", occ.pth))
    sp.df <- data.frame("species" = species)
    sp.df[, c("latitude_mean", "longitude_mean")] <- rep(NA, length(sp.df$species))
    # Extracting ...
    for(i in 1:length(sp.df$species)){
        print(paste0(i, ": ", sp.df$species[i]))
        occ.file <- paste0(occ.pth, str_replace(sp.df$species[i], " ", "_"), ".csv")
        if(file.exists(occ.file)){
          occ.df <- read.csv(occ.file, header = T)
          coord <- coordinates(occ.df[, c("longitude", "latitude")])
          meanValues <- apply(coord, 2, FUN = mean, na.rm = T)
            sp.df[i, "longitude_mean"] <- meanValues[1]
            sp.df[i, "latitude_mean"] <- meanValues[2]
          }
        }
  return(sp.df)
}

## Mean longlat calculation and merge with ht

ht <- read.csv("data/HydraTRY_Pablo_2020.csv", header = T)  # Includes three principal components

longLat <- meanLongLatFun(species = ht$Species, occ.pth = "data/species_data/occurrences/processed_data_outl/")
longLat$Species <- longLat$species
ht_longlat <- join(ht, longLat, by ="Species")
head(ht_longlat)
length(ht_longlat$Species)

# Save

write.csv(ht_longlat, "data/HydraTRY_Pablo_2020.csv", row.names = F)
print("====> data/HydraTRY_Pablo_2019.csv")
