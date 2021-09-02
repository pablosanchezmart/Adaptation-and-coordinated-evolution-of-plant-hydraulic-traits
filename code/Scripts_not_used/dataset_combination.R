#################################### DATASET MERGING ##########################################################

# Sanchez Pablo
# 25-02-2018

### PACKAGES ####
# install.packages("AgrgregateR", dependencies = T)
library(dplyr)

### DATAFRAMES AND PARAMETERS ####

# New Species dataframe (newSp.fun)
newSpp.df <- data.frame(1)

# New values dataframe 
newValues.df <- data.frame(1)

# To aggregate values
newSppAgr.df <- data.frame(1)

# Counter to add new values to newValues.df
n <- 0

## To aggregate dataframes
#Numerical variables
numVar <- c("Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "Pitlp", "AlAs", "WD", "Vdia", "VD", "VLmax", "ALT", "AI")
# Categorical variables
catVar <- c("Life.form", "Leaf.form", "Soiltype", "biome")


### FUNCTIONS ####

## New species function
#To check if there is any new specie for our database.
newSp.fun <- function(mainDB, database2){
  for(i in 1:length(database2$Species)){                                                             # Bucle for all the recorded
    sp <- database2$Species[i]                                                                       # Select a specie of the exogenous dataset
    # New specie?
    if(length(mainDB[which(as.character(mainDB$Species) == as.character(sp)), ]$Species) == 0){      # Check if the specie is new                                                                                                         # New specie row
      levels(newSpp.df) <- as.factor(sp)                                                                        # Add specie as a factor level
      newSpp.df <- bind_rows(newSpp.df, database2[which(database2$Species ==sp), ])                      # Bind new row to the new species data frame
      cat("New specie:", sp, "\n")
    } else {
      next()
    }
  } # end of bucle
  newSpp.df$Species <- factor(newSpp.df$Species)
  return(newSpp.df)
} # end of function

## Filling gaps function
#To fill gaps in the existing database variables
fillDB.fun <- function(mainDB, database2, var) {
  for(i in 1:length(database2$Species)){
    sp <- mainDB$Species[i]  
      if(is.na(mainDB[which(mainDB$Species == sp), var]) && length(database2[which(database2$Species == sp), var]) != 0){                 # HydraTRY: NA, new.Hmax: no NA (Hmax value)
        mainDB[which(mainDB$Species == sp), var] <- database2[which(database2$Species == sp), var]                                        # Assign new Hmax value to HtdraTRY
        cat("New record for: ", var, " SPECIES: ", sp, "\n")                                                                              # Print that we have a new record
      } else {                                                                                                                            # If HydraTRY: value and/or new.Hmax: NA
        if(!is.na(mainDB[which(mainDB$Species == sp), var]) && length(database2[which(database2$Species == sp), var]) != 0){              # If we have a new value for a specie that already had in HydraTRY
          if(!is.na(database2[which(database2$Species == sp), var])){                                                                      # To avoid NAs problems
          if(mainDB[which(mainDB$Species == sp), var] != database2[which(database2$Species == sp), var]) {                                # If values are different
            cat("New measure for", sp, ".\n", var, "Main database record :", mainDB[which(mainDB$Species == sp), var], 
                   "\n New value :", database2[which(database2$Species == sp), var], "\n")                               
            n <- n + 1                                                                                                                    # Counter for the record of the new measures
            newValues.df[n, "Species"] <- sp
            newValues.df[n, var] <- database2[which(database2$Species == sp), var]
          } # Close of new measure
        } # Close NA's avoiding
      }  else {                                                                                                                          # If both values are the same, next.
          next()                                                      
      } # End repeated measure
  } # end of else
} # End of bucle
return(newValues.df)
return(mainDB)
} # End of function

## Mode function (for categorical variables)
mode.fun <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

###################### CODE TO COMPLETE DATA FRAME ####

###################### LIU ET AL. 2019 ################
setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")

### NEW SPECIES AND VARIABLES ####

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
length(levels(liuData$Species))
length(liuData$Species)
summary(liuData)

# Variables to refill (name/description: Liu variable name (units) = HydraTRY variable name (units))
# *New to new variables interesting to incorpore.

# Maximum heightHmax = Hmax
# Actual height: Hact (m) = *New
# Xylem conductivity: Ks (Kg m-1 MPa-1 s-1) = Ks (Kg m-1 MPa-1 s-1)
# Leaf hydraulic conductivity: Kl (10^-4 Kg m-1 MPa-1 s-1) = Kl (Kg m-1 MPa-1 s-1)
# Xylem tension at which 50% of the maximum hydraulic conductivity is lost: P50 (MPa) = P50 (MPa)
# Minimum predawn xilem preassure potencial: Ppd (MPa) = MinWP_pd (MPa)
# Minimum midday predawn xilem preassure potencial: Pmd (MPa) = MinWP_md
# Leaf turgor loss point: TLP (MPa) = Pitlp (MPa)
# Leaf area to sapwood area ratio AlAs (m2 cm-2) = 1/Hv (m2 cm-2)
# Sapwood density/ Wood density: WD (g cm-3) = WD (wood density) (g cm-3) 

# Aridity Index (MAP/PET, http://www.cgiar-csi.org): AI = *New

### Preparation of database (same names)

names(liuData)[names(liuData) == "Ppd"] <- "MinWP_pd"
names(liuData)[names(liuData) == "Pmd"] <- "MinWP_md"
names(liuData)[names(liuData) == "TLP"] <- "Pitlp"
names(liuData)[names(liuData) == "TLP"] <- "Pitlp"
names(liuData)[names(liuData) == "Life form"] <- "Life.form"
names(liuData)[names(liuData) == "Leaf form"] <- "Leaf.form"
names(liuData)[names(liuData) == "Biome"] <- "biome"
names(liuData)[names(liuData) == "Group"] <- "group"
names(liuData)[names(liuData) == "Family"] <- "family"
names(liuData)[names(liuData) == "Genus"] <- "genus"

# Note that Kl needs to be divided by 10^4 to have same units (Liu_etal-2019_Data.csv already changed)
liuData$Kl <- liuData$Kl/(10^4)

# Change level names for Biome factor
levels(ht$biome)
levels(liuData$biome)
liuData$biome <- plyr::revalue(liuData$biome, c("BOR"= "Boreal forest", "DES" = "Desert", "TMR" = "Temperate rain forest", "TMS" = "Temperate forest", "TRR" = "Tropical rain forest",
       "TRS" = "Tropical seasonal rainforest","WDS" = "Woodland/shrubland"))

### Add empty new variables to HydraTRY
ht$Hact <- numeric(length(ht$Species))
ht$Hact <- NA
ht$AlAs <- numeric(length(ht$Species))
ht$AlAs <- NA
ht$Vdia <- numeric(length(ht$Species))
ht$Vdia <- NA
ht$VD <- numeric(length(ht$Species))
ht$VD <- NA
ht$VLmax <- numeric(length(ht$Species))
ht$VLmax <- NA
ht$ALT <- numeric(length(ht$Species))
ht$ALT <- NA
ht$Life.form <- factor(length(ht$Species))
ht$Life.form <- NA
ht$Leaf.form <- factor(length(ht$Species))
ht$Leaf.form <- NA
ht$Soiltype <- factor(length(ht$Species))
ht$Soiltype <- NA

### New species of Liu et al. data compared to HydraTRY
newSpp.df <- newSp.fun(ht, liuData)


# Aggregate data
SppNumV.df <- aggregate(newSpp.df[, numVar], by = list(newSpp.df$Species), FUN = mean, na.action = na.omit)
SppCatV.df <- aggregate(newSpp.df[, catVar], by = list(newSpp.df$Species), FUN = mode.fun)

# Merge categorical and numerical variables
newSppAgr.df <- merge(SppNumV.df, SppCatV.df, by = "Group.1")
names(newSppAgr.df)[names(newSppAgr.df) == "Group.1"] <- "Species"
length(newSppAgr.df$Species)
# 577 new species

### Hv calculation ###
# For the new spp to incorporate
newSppAgr.df$Hv <- 1/newSppAgr.df$AlAs
# Also for Liu
liuData$Hv <- 1/liuData$AlAs
### Add species to HydraTRY ###

ht_liu <- bind_rows(newSppAgr.df, ht)

ht_liu$X <- NULL
ht_liu$Species <- as.factor(ht_liu$Species)
names(ht_liu)
summary(ht_liu)
summary(liuData)
length(ht_liu$Species)
### Save new variables in HydraTRY dataset ###

# write.csv(ht_liu, "HydraTRY_Pablo_logData_HmaxCompl_Liu.csv", row.names = FALSE)

### NEW VALUES FOR EXISTING VARIABLES ####

### Aggregate Liu dataset

liuDataNum <- aggregate(liuData[, numVar], by = list(liuData$Species), FUN = mean, na.action = na.omit)
liuDataCat <- aggregate(liuData[, catVar], by = list(liuData$Species), FUN = mode.fun)

# Merge categorical and numerical variables
liuData <- merge(liuDataNum, liuDataCat, by = "Group.1")

names(liuData)[names(liuData) == "Group.1"] <- "Species"
length(liuData$Species)

# Change species factors to characters
ht_liu$Species <- as.character(ht_liu$Species)  # something is wrong when working with factors, that is a easy solution.
liuData$Species <- as.character(liuData$Species)

# Variables to extract
variables <- names(ht_liu[, -1])
# Bucle of extraction
for(variable in variables){
  cat("Scrapping", variable, "\n")
  fillDB.fun(ht_liu, liuData, variable)
}

# Checking
ht_liu$Species <- as.factor(ht_liu$Species)
ht_liu$Life.form <- as.factor(ht_liu$Life.form)
ht_liu$Leaf.form <- as.factor(ht_liu$Leaf.form)
ht_liu$Soiltype <- as.factor(ht_liu$Soiltype)
ht_liu$biome <- as.factor(ht_liu$biome)


head(ht_liu)
length(ht_liu$Species)
summary(ht_liu)
ht_liu <- arrange(ht_liu, Species)  # Order alphabetically
head(ht_liu)
levels(ht_liu$biome)  # Maybe we can merge some of this biomes
names(ht_liu)
length(ht_liu)
ht_liu <- ht_liu[,c("Species", "genus", "family", "order", "group", "Growth.form", "Life.form", "Leaf.form", "biome", "Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "Pitlp", "AlAs", "WD", "Vdia", "VD", "VLmax", "ALT", "AI",
               "Hv", "HSM", "PItlp", "WUE", "SLA", "Aarea", "Amass", "LL", "pLL", "N", "Narea", "Rd", "lat_mean", "long_mean", "MAT_mean", "MAP_mean", "tmin", "tmax", "srad", "wind", "vapr", "TMean.year", "Diurnal.T.Range", 
               "Isothermality", "Seasonality", "TMax.warmest", "TMin.Coldest" , "TRange", "TMean.Wettest", "TMean.Driest", "TMean.Warmest", "TMean.Coldest", "Prec.year", "Prec.Wettest", "Prec.Driest", "Prec.Seasonality", "Prec.WetQ",
               "Prec.DryQ", "Prec.WarmQ", "Prec.ColdQ", "PPET", "AbsoluteDepthBedRock_1km_ll", "PredictedClassTAXOUSDA_1km_ll",  "PredictedClassTAXNWRB_1km_ll", "Soiltype")]
head(ht_liu)

#write.csv(ht_liu, "HydraTRY_Pablo_logData_HmaxCompl_Liu.csv", row.names = F)

newValues.df  # If we want to see new values, run HydraTRY_Pablo_logData_HmaxCompl.csv again from above.
