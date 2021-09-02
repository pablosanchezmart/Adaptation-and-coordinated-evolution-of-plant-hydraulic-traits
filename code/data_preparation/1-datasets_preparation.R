##### PREPARE DATASETS TO JOIN HT ###############################################################################

# P Sanchez-Martinez

print("Preparing datasets to join Ht")

#### PACKAGES --------------------------------------------------------------------------------------------- ####
library(stringr)
library(knitr)
library(dplyr)
library(plyr)
#vpd 
library(humidity) # for saturation vapor preassure calculation from T
#https://cran.r-project.org/web/packages/humidity/humidity.pdf

#### FUNCTIONS ------------------------------------------------------------------------------------------- ####

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

#### Hmax --------------------------------------------------------------------------------------------- ####

hmax <- read.csv("data/Others/HydraTRY_Pablo_HmaxCompl.csv", header = T)
hmax <- hmax[, c("Species", "Hmax") ]
write.csv(hmax, "data/processed_datasets/hmax.csv", row.names = F)

#### LIU ET AL. 2019 ----------------------------------------------------------------------------------- ####

## Charge Liu et al. 2019 database...
liuData <- readxl::read_excel("data/Others/Liu_etal-2019_Data.xlsx", sheet = "full_data")
length(liuData$Species)

# Select variables
names(liuData)
liuData <- liuData[, -c(1, 2, 9, 10, 11, 25, 30, 31, 33)]
liuData$Species <- as.factor(liuData$Species)
names(liuData)

### Preparation of database (same names)
names(liuData)[names(liuData) == "Ppd"] <- "MinWP_pd"
names(liuData)[names(liuData) == "Pmd"] <- "MinWP_md"
names(liuData)[names(liuData) == "TLP"] <- "PItlp"
names(liuData)[names(liuData) == "Life form"] <- "life.form"
names(liuData)[names(liuData) == "Leaf form"] <- "leaf.form"
names(liuData)[names(liuData) == "Growth.form"] <- "growth.form"
names(liuData)[names(liuData) == "Biome"] <- "biome"
names(liuData)[names(liuData) == "Group"] <- "group"
names(liuData)[names(liuData) == "Family"] <- "family"
names(liuData)[names(liuData) == "Genus"] <- "genus"

# Note that Kl needs to be divided by 10^4 to have same units (Liu_etal-2019_Data.csv already changed)
summary(liuData$Kl)
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
numVar <- c("Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "PItlp", "AlAs", "WD", "Vdia", "VD", "VLmax", "Hv", "AI", "MAP", "MAT")
# Categorical variables
catVar <- c("life.form", "leaf.form", "biome", "family", "genus", "group")

liuDataNum <- aggregate(liuData[, numVar], by = list(liuData$Species), FUN = mean.fun)  # Had problems with the previus mean fun (not account for NAs omission)
liuDataCat <- aggregate(liuData[, catVar], by = list(liuData$Species), FUN = mode.fun)

# Merge categorical and numerical variables
liuData <- merge(liuDataNum, liuDataCat, by = "Group.1")

names(liuData)[names(liuData) == "Group.1"] <- "Species"
names(liuData)

# Save prepared dataset
write.csv(liuData, "data/processed_datasets/liu_etal-2019_DataJoin.csv", row.names = F)

#### GS DATASET ---------------------------------------------------------------------------------------- ####

gs <- readxl::read_excel("data/Others/gsDataset.xlsx")
gs$gmax <- as.numeric(gs$`gmax (mol m-2 s-1)`)
gs$Species <- str_replace_all(gs$Species, "_", " ")
gs$Species <- word(gs$Species, 1, 2)
gs$Species <- as.factor(gs$Species)
length(unique(gs$Species))

## Agregate data
gsData <- aggregate(gs[, "gmax"], by = list(gs$Species), FUN = mean.fun)
names(gsData)[names(gsData) == "Group.1"] <- "Species"

# Save prepared dataset
write.csv(gsData, "data/processed_datasets/gsDataset.csv", row.names = F)

#### HV dataset ------------------------------------------------------------------------------------------ ####

## Charge HV dataset
hvData <- read.csv("data/Others/fin_HV_Pablo.csv", header = T)
head(hvData)
names(hvData)[names(hvData) == "Decid_Ever"] <- "leaf.form"
summary(hvData)

write.csv(hvData, "data/processed_datasets/life_form.csv", row.names = F)