##### ENVIRONMENTAL VARIABLES EXTRACTION #############################################################################

# P Sanchez-Martinez

print("Extracting environmental variables and including them in dataset...")

### PACKAGES --------------------------------------------------------------------------------------------####
# install.packages("devtools")
# devtools::install_github("ropensci/CoordinateCleaner")

library(sp)
library(stringr)
library(ALA4R)
library(rgbif)
library(raster)
library(futile.logger)
library(doBy)
library(knitr)
library(dplyr)
library(dismo)
library(humidity)

library(CoordinateCleaner)
library(countrycode)

#### FUNCTIONS ----------------------------------------------------------------------------------------------- ####

# Annualize variables
anuMeanWCFun <- function(vars =  c("tmin", "tmax", "srad", "wind", "vapr"), res = "30s"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/env_data/worldclim/wc", res, "/", var, "_mean_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/env_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = mean, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      removeTmpFiles(h=1) # remove temporal raster files older than 1 hour
    }
  }
}

anuMaxWCFun <- function(vars =  c("wind"), res = "30s"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/env_data/worldclim/", var, "_max_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/env_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = max, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      removeTmpFiles(h=1) # remove temporal raster files older than 1 hour
    }
  }
}

anuMinWCFun <- function(vars =  c("srad"), res = "30s"){
  for(var in vars){
    print(var)
    file.name <- paste0("data/species_data/env_data/worldclim/", var, "_min_annual.tif")
    if(file.exists(file.name)){
      next()
    } else {
      var.list <- list.files(path = paste0("data/species_data/env_data/worldclim/wc", res, "/"), pattern = paste0("wc2.1_", res, "_", var), full.names = T)
      var.st <- stack(var.list)
      anuVar.r <- raster::clusterR(x = var.st, overlay, args = list(fun = min, na.rm = T))
      names(anuVar.r) <- var
      writeRaster(anuVar.r, file.name, overwrite = T)
      removeTmpFiles(h=1) # remove temporal raster files older than 1 hour
    }
  }
}

# Occurrences download

gbifOccFun <- function(species){
  for(sp in species){
    print(sp)
    output.file <- paste0("data/species_data/occurrences/raw_data/", str_replace(sp, " ", "_"), ".csv")
    if(file.exists(output.file)){
      next()
    } else {
      # ocurences download
      gbif.data <- occ_search(scientificName = sp)                                                  # That one takes time (to 20 min /sp)
      if(gbif.data$meta$count == 0){
        next()
      }
      if(!exists("decimalLatitude", where = gbif.data$data)){
        next()
      }
      if(all(is.na(gbif.data$data$decimalLatitude))){
        next()
      } else {
        gbif.df <- na.omit((as.data.frame(gbif.data$data[, c("decimalLatitude", "decimalLongitude")])))
        occ.df <- data.frame("species" = sp,"longitude" = gbif.df$decimalLongitude, 
                             "latitude" = gbif.df$decimalLatitude)
        write.csv(occ.df,  output.file, row.names = F)
        print(paste0("Occurrences saved: ", output.file))
      }
    }
  }
}

ALAOccFun <- function(species){
  for(sp in species){
    print(sp)
    output.file <- paste0("data/species_data/occurrences/raw_data/", str_replace(sp, " ", "_"), ".csv")
    if(file.exists(output.file)){
      next()
    } else {
      # ocurences download
      occ.df <- ALA4R::occurrences(sp, email = "p.sanchez@creaf.uab.cat", download_reason_id = 7)                                                  # That one takes time (to 20 min /sp)
      if(all(is.na(occ.df$longitude))){
        next()
      } else {
        write.csv(occ.df, output.file, row.names = F)
        print(paste0("Occurrences saved: ", output.file)) 
      }
    }
  }
}

cleanOccFun <- function(species){
  i <- 1
  NASp <- c()
  # Filter results
  filter.file <- "data/species_data/occurrences/processed_data/filtering_results.csv"
  if(file.exists(filter.file)){
    filter.df <- read.csv(filter.file, header = T)
    filter.df$species <- as.character(filter.df$species)
    n <- length(filter.df[, 1]) + 1
  } else {
    filter.df <- data.frame("species" = character(),"raw" = numeric(), "processed" = numeric(), "filtered" = numeric())
    filter.df$species <- as.character(filter.df$species)
    n <- 1
  }
  # Occurrences cleaning
  for(sp in species){
    print(sp)
    input.file <-  paste0("data/species_data/occurrences/raw_data/", str_replace(sp, " ", "_"), ".csv")
    output.file <- paste0("data/species_data/occurrences/processed_data/", str_replace(sp, " ", "_"), ".csv")
    # File exists?
    if(!file.exists(output.file)){
      if(file.exists(input.file)){
        occ.df <- read.csv(input.file, header = T)
        # Data exists?
        if(all(is.na(occ.df$longitude))){
          file.remove(input.file)
          cat("deleting empty file: ", input.file, "\n")
          NASp[i] <- sp
          i <- i + 1
        } else {
          # Occurrences cleaning
          clean <- clean_coordinates(x = occ.df, lon = "longitude", lat = "latitude",
                                       species = "species",
                                       tests = c("capitals", "duplicates","equal","gbif", "institutions", "seas", "zeros"), value = c("clean"), verbose = F) # most test are on by default
          write.csv(clean, output.file, row.names = F)
          filter.df[n, "species"] <- as.character(sp)
          filter.df[n, "raw"] <- length(occ.df$species)
          filter.df[n, "processed"] <- length(clean$species)
          filter.df[n, "filtered"] <- length(occ.df$species) - length(clean$species)
          n <- n + 1
        }
      } else {
        NASp[i] <- as.character(sp)
        i <- i + 1
      }
    }
  }
  write.csv(filter.df, filter.file, row.names = F)
  rslts <- list()
  rslts$NASp <- NASp
  rslts$filter.df  <- filter.df
  return(rslts)
}

cleanOutlFun <- function(species){
  i <- 1
  NASp <- c()
  # Filter results
  filter.file <- "data/species_data/occurrences/processed_data/filtering_results.csv"
  if(file.exists(filter.file)){
    filter.df <- read.csv(filter.file, header = T)
    filter.df$species <- as.character(filter.df$species)
    n <- length(filter.df[, 1]) + 1
  } else {
    filter.df <- data.frame("species" = character(),"raw" = numeric(), "processed" = numeric(), "filtered" = numeric())
    filter.df$species <- as.character(filter.df$species)
    n <- 1
  }
  # Occurrences cleaning
  for(sp in species){
    print(sp)
    input.file <-  paste0("data/species_data/occurrences/processed_data/", str_replace(sp, " ", "_"), ".csv")
    output.file <- paste0("data/species_data/occurrences/processed_data_outl/", str_replace(sp, " ", "_"), ".csv")
    # File exists?
    if(!file.exists(output.file)){
      if(file.exists(input.file)){
        occ.df <- read.csv(input.file, header = T)
        # Data exists?
        if(all(is.na(occ.df$longitude))){
          file.remove(input.file)
          cat("deleting empty file: ", input.file, "\n")
          NASp[i] <- sp
          i <- i + 1
        } else {
          # Occurrences cleaning
          clean <- cc_outl(occ.df, lon = "longitude", lat = "latitude", species = "species", method = "quantile", value = "clean")
          write.csv(clean, output.file, row.names = F)
          filter.df[n, "species"] <- as.character(sp)
          filter.df[n, "raw"] <- length(occ.df$species)
          filter.df[n, "processed"] <- length(clean$species)
          filter.df[n, "filtered"] <- length(occ.df$species) - length(clean$species)
          n <- n + 1
        }
      } else {
        NASp[i] <- as.character(sp)
        i <- i + 1
      }
    }
  }
  write.csv(filter.df, filter.file, row.names = F)
  rslts <- list()
  rslts$NASp <- NASp
  rslts$filter.df  <- filter.df
  return(rslts)
}

# Extract environmental variables species mean values
extrSpFun <- function(species = ht$Species, pth = paste0("data/species_data/soil_data/"), occ.pth = "data/species_data/occurrences/processed_data_outl/", vars = c("sg"), sttstc = mean, prev.rslts = F){
  # Variables
  for(var in vars){
    print(paste0("extracting data from ", pth))
    var.list <- list.files(path = pth, pattern = ".tif", full.names = T)
    print(var.list)
    st <- stack(var.list)
    varNames <- names(st)
    # Load previous results, if exist
    if(isTRUE(prev.rslts)){
      sp.df <- read.csv(paste0("data/species_data/ht_species_env_", var, ".csv"), header = T)
    } else {
      sp.df <- data.frame("species" = species)
      sp.df[, varNames] <- rep(NA, length(sp.df$species))
    }
    # Extracting ...
    for(i in 1:length(sp.df$species)){
      if(any(is.na(sp.df[i, varNames]))){
        print(paste0(i, ": ", sp.df$species[i]))
        occ.file <- paste0(occ.pth, str_replace(sp.df$species[i], " ", "_"), ".csv")
        if(file.exists(occ.file)){
          occ.df <- read.csv(occ.file, header = T)
          coord <- coordinates(occ.df[, c("longitude", "latitude")])
          extr.data <- raster::extract(x = st, y = coord)
          varValues <- apply(extr.data, 2, FUN = sttstc, na.rm = T)
          for(varName in varNames){
            sp.df[i, varName] <- varValues[varName] 
          }
          write.csv(sp.df, paste0("data/species_data/ht_species_env_", var, ".csv"), row.names = F)
          removeTmpFiles(h=0.5) # remove temporal raster files older than 0.5 hour
        } else {
          write.csv(sp.df, paste0("data/species_data/ht_species_env_", var, ".csv"), row.names = F)
        }
      } else {
        next()
      }
    }
  }
  return(sp.df)
}

# To divide temperature variables per 10
tempVarCorrection.fun <- function(df, vars){
  for(i in vars){
    df[[i]] <- df[[i]]/10
  } # end bucle
  return(df)
} # end function

## VPDmax warmest month calculation
vpdMaxFun <- function(df = env, tmaxVars = c(23:34)){
  df$VPDMax <- numeric(length(df[,1]))
  print("Calculating species VPDmax warmest month ...")
  for(i in 1:length(df$species)){
    if(all(is.na(df[i, tmaxVars]))){
      df$VPDMax[i] <- NA
    } else {
      tmax.month <- which.max(df[i , tmaxVars])                                                     # Warmest month
      if(tmax.month < 10){
        tmax.month <- paste0("0", as.character(tmax.month))
      }
      vapr.month <- paste0("vapr_", tmax.month)                                                  # vapr (kPa) of the warmest month           
      tmax.month <- paste0("tmax_", tmax.month)                                                  # Tmax value warmest month
      actVapr <- df[i, vapr.month]                                                               # vapr value warmest month
      satVp <- SVP(df[i, tmax.month], isK = F, formula = c("Clausius-Clapeyron", "Murray"))      # SatVP (hPa) (Tmax warmest month)
      df$VPDMax[i] <- (satVp/10) - actVapr                                                       # VPDmax = vapr(warmest month) - Tmax (warmest month)
    }
  }
  print("Removing monthly tmax and vapr variables...")
  df <- df[, -tmaxVars]
  vaprVars <- paste0("vapr_0", 1:9)
  vaprVars <- c(vaprVars, paste0("vapr_", 10:12))
  vaprVars <- which(names(df) %in% vaprVars)
  df <- df[, -vaprVars]
  return(df)
}# end of function

### Database -------------------------------------------------------------------------------------------------------------- ####

ht <- read.csv("data/ht/HydraTRY_Pablo_traits_taxscrub.csv", header = T)  # Last version

#### PREPARE ENVIRONMENTAL DATA ------------------------------------------------------------------------------------------- ####
print("preparing environmental data...")

# Annualize WC climate data

beginCluster()
anuMaxWCFun(vars = c("wind"), res = "30s")
anuMinWCFun(vars = c("srad"), res = "30s")
anuMeanWCFun(vars = c("tmin"), res = "30s")
endCluster()


#### OCCURRENCES --------------------------------------------------------------------------------------------- ####

gbifOccFun(species = ht$Species)
ALAOccFun(species = ht$Species)

## Cleaning
ht1 <- ht[1:1047, ]
clean.rslts <- cleanOccFun(species = ht1$Species)
clean.rslts$filter.df

outl.rslts <- cleanOutlFun(species = ht1$Species)

### DATA EXTRACTION ------------------------------------------------------------------------------------------ ####

beginCluster()
ht_sp.df <- extrSpFun(species = ht$Species, pth = paste0("data/species_data/env_data/soil_data/new/"), occ.pth = "data/species_data/occurrences/processed_data_outl/", vars = c("new_sg"),
                      sttstc = mean, prev.rslts = T)
endCluster()

r <- read.csv("data/species_data/ht_species_env_wc.csv", header = T)
names(r)

#### CLIMATE DATASET ----------------------------------------------------------------------------------- ####

wc <- read.csv("data/species_data/ht_species_env_wc.csv", header = T)
sg <- read.csv("data/species_data/ht_species_env_sg.csv", header = T)
ai <- read.csv("data/species_data/ht_species_env_ai.csv", header = T)

env <- merge(wc, sg, by = "species")
env <- merge(env, ai, by = "species")
names(env)
envNames <- str_remove_all(names(env), "wc2.1_30s_")
names(env) <- envNames

#### VPDmax calculation -------------------------------------------------------------------------------- ####

env.vpd <- vpdMaxFun(df = env, tmaxVars = c(23:34))
names(env.vpd)

env.vpd <- dplyr::rename(env.vpd, 
                         # Worldclim
                         Species = species, MAT = bio_1, Diurnal_TRange = bio_2, Isotermality = bio_3, TSeasonality = bio_4, TMax.Warmest = bio_5,
                         TMin.Coldest = bio_6, TRange = bio_7, TMean_Wettest = bio_8, TMean_Driest = bio_9, TMean_Warmest = bio_10, TMean_Coldest = bio_11,
                         AP = bio_12, Prec_Wettest = bio_13, Prec_Driest = bio_14, Prec_Seasonality = bio_15, Prec_WetQ = bio_16, Prec_DryQ = bio_17, Prec_WarmQ = bio_18,
                         Prec_ColdQ = bio_19, windMax = wind_max_annual, wind = wind_mean_annual, vapr = vapr_mean_annual, srad = srad_mean_annual, sradMin = srad_min_annual, 
                         tmax = tmax_mean_annual, tmin = tmin_mean_annual,
                         # Soilgrids
                         ABDRock = AbsoluteDepthBedRock_1km_ll, clay_60cm = ClayContent_percentage_60cm_1km, ph_60cm = PH_HOX_60cm_1kmll,
                         silt_60cm = SiltContent_percentage_60cm_1km, sand_60cm = SandContent_percentage_60cm_1km, SWC_200cm = SaturatedWaterContent_volumetricFraction_200cm_1km,
                         # Aridity index (Moisture index)
                         MI = AI)

names(env.vpd)
summary(env.vpd)


#### VARIABLE CALCULATION/RECALCULATION ---------------------------------------------------------------- ####

# For temperature and PH vars we need to divide the values by 10
divVars <- c("ph_60cm", "MI")

env.vpd <- tempVarCorrection.fun(env.vpd, divVars)
summary(env.vpd)

# Save prepared dataset
write.csv(env.vpd, "data/processed_datasets/environmental_data.csv", row.names = F)

### MERGE HYDRATRY WITH CLIM DATASET ####

names(ht)
names(env.vpd)

## Combine datasets
ht_env <- merge(ht, env.vpd, by = "Species", all = TRUE)
names(ht)
names(ht_env)
summary(ht_env)

# Checking
names(ht_env)
length(which(!ht_env$Species %in% ht$Species))
ht_env[which(!ht_env$Species %in% ht$Species), "Species"]
length(ht_env$Species)
length(ht$Species)
summary(ht_env)

ht_env_fin <- ht_env %>% filter(ht_env$Species %in% ht$Species)
length(ht_env_fin$Species)
summary(ht_env_fin)

# Save complete dataframe

write.csv(ht_env_fin, "data/HydraTRY_Pablo_2020.csv", row.names = F)
