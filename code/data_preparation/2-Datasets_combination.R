##### DATA COMBINATION ####################################################################################

# P Sanchez-Martinez

print("Combining traits datasets...")

### PACKAGES ----------------------------------------------------------------------------------------- ####
library(plyr)

# Original dataset
ht <- read.csv("data/ht/HydraTRY_Pablo_original.csv", header = T)

### FUNCTIONS ----------------------------------------------------------------------------------------- ####

## Combination function
comb.fun <- function(df, var){
  var.x <- paste(var, ".x", sep = "")
  var.y <- paste(var, ".y", sep = "")
  # choose default value from 'x' in the case that both '.x' and '.y' have values
  df[!is.na(df[, var.x]) & !is.na(df[, var.y]), var] <- df[!is.na(df[, var.x]) & !is.na(df[, var.y]), var.x]
  # replace NA values in 'x' with values from 'y' and vice-versa
  df[is.na(df[, var.x]), var] <- df[is.na(df[, var.x]), var.y]
  df[is.na(df[, var.y]), var] <- df[is.na(df[, var.y]), var.x]
  # Delete ".x" and ".y" columns
  df[, var.x] <- NULL
  df[, var.y] <- NULL
  return(df)
}

invComb.fun <- function(df, var){
  var.x <- paste(var, ".x", sep = "")
  var.y <- paste(var, ".y", sep = "")
  # choose default value from 'y' in the case that both '.x' and '.y' have values
  df[!is.na(df[, var.x]) & !is.na(df[, var.y]), var] <- df[!is.na(df[, var.x]) & !is.na(df[, var.y]), var.y]
  # replace NA values in 'x' with values from 'y' and vice-versa
  df[is.na(df[, var.x]), var] <- df[is.na(df[, var.x]), var.y]
  df[is.na(df[, var.y]), var] <- df[is.na(df[, var.y]), var.x]
  # Delete ".x" and ".y" columns
  df[, var.x] <- NULL
  df[, var.y] <- NULL
  return(df)
}

#### HMAX --------------------------------------------------------------------------------------------- ####

# Data
hmax <- read.csv("data/processed_datasets/hmax.csv", header = T)

# Combine datasets

ht_hmax <- merge(ht, hmax, by = "Species", all = T)
ht_hmax <- comb.fun(ht_hmax, "Hmax")

# Checking
with(ht, length(which(!is.na(Hmax)))) # ht data
with(hmax, length(which(!is.na(Hmax)))) # hmax data
with(ht_hmax, length(which(!is.na(Hmax)))) # Comb dataset data
length(ht_hmax$Species)
names(ht_hmax)

#### LIU ETAL 2019 ------------------------------------------------------------------------------------- ####

## Data
liuData <- read.csv("data/processed_datasets/liu_etal-2019_DataJoin.csv", header = T)

## Combine datasets

# First join new colunmns
ht_hmax_liu <- join(ht_hmax, liuData, by = "Species", type = "full")
names(ht_hmax)
names(ht_hmax_liu)

# Combine datasets
ht_hmax_liu <- merge(ht_hmax_liu, liuData, by = "Species", all = T, sort = T)
names(ht_hmax_liu)
variables <- names(liuData) # Combine only shared variables
variables <- variables[-1]  # Drop Species

for(variable in variables){
  ht_hmax_liu <- comb.fun(ht_hmax_liu, variable)
}

# Arrange alphabetically
ht_hmax_liu <- arrange(ht_hmax_liu, Species)

## Checking
names(ht)
names(ht_hmax_liu)
length(ht_hmax_liu$Species)
summary(ht_hmax_liu$Kl)

# Non shared specie
ht_hmax_liu[which(ht_hmax_liu$Species == "Zanthoxylum avicennae"), ]
ht[which(ht$Species == "Zanthoxylum avicennae"), ]
liuData[which(liuData$Species == "Zanthoxylum avicennae"), ]
# Shared specie with complete values (e.g. Hmax)
ht_hmax_liu[which(ht_hmax_liu$Species == "Acacia havilandiorum"), ]
ht[which(ht$Species == "Acacia havilandiorum"), ]
liuData[which(liuData$Species == "Acacia havilandiorum"), ]
# Shared specie fill gap (e.g. Ks)
ht_hmax_liu[which(ht_hmax_liu$Species == "Acacia greggii"), ]
ht[which(ht$Species == "Acacia greggii"), ]
liuData[which(liuData$Species == "Acacia greggii"), ]

with(ht_hmax_liu, length(which(!is.na(Hmax))))
with(ht_hmax_liu, length(which(!is.na(Hv))))
with(ht_hmax_liu, length(which(!is.na(Ks))))
length(ht_hmax_liu$Species)

summary(ht_hmax_liu)

### MERGE HYDRATRY WITH GS DATASET --------------------------------------------------------------------------------- ####

gs <- read.csv("data/processed_datasets/gsDataset.csv", header = T)
names(gs)
## Combine datasets
ht_hmax_liu_gs <- merge(ht_hmax_liu, gs, by = "Species", all = TRUE)

# Checking
names(ht_hmax_liu_gs)
with(ht_hmax_liu_gs, length(which(!is.na(gmax))))
length(ht_hmax_liu_gs$Species)
length(ht_hmax_liu_gs)
summary(ht_hmax_liu_gs)

### MERGE HYDRATRY WITH DECIDUOUS/EVERGREEN DATASET ####

# Data
deData <- read.csv("data/processed_datasets/life_form.csv", header = T)

# Combine datasets
summary(ht_hmax_liu_gs$leaf.form)
summary(deData$leaf.form)

ht_hmax_liu_gs_de <- merge(ht_hmax_liu_gs, deData, by = "Species", all.x = T)
names(ht_hmax_liu_gs_de)
length(ht_hmax_liu_gs_de$Species)

ht_hmax_liu_gs_de <- invComb.fun(ht_hmax_liu_gs_de, "leaf.form")

# Checking
summary(ht_hmax_liu_gs_de$leaf.form)
summary(ht_hmax_liu_gs_de)
length(ht_hmax_liu_gs_de$Species)

# Write combined Dataset
write.csv(ht_hmax_liu_gs_de, "data/ht/HydraTRY_Pablo_traits.csv", row.names = F)

### DATASET FILLING USING KNOWN RELATIONSHIPS BETWEEN VARIABLES ####

ht <- read.csv("data/ht/HydraTRY_Pablo_traits.csv", header = T)

## Hv-Ks-Kl relationship
summary(ht$Hv) # 1108 NAs
summary(ht$Ks) # 1238 NAs
summary(ht$Kl) # 1478 NAs

# Kl = Hv * Ks
names(ht)[names(ht) == "Kl"] <- "Kl.x"
ht$Kl.y <- (ht$Hv / 10000) * ht$Ks
summary(ht$Kl.x)
summary(ht$Kl.y)
ht <- comb.fun(ht, "Kl")
summary(ht$Kl)  # 1396 NAs

# Hv = Kl/Ks
names(ht)[names(ht) == "Hv"] <- "Hv.x"
ht$Hv.y <- (ht$Kl / ht$Ks) * 10000
summary(ht$Hv.x)
summary(ht$Hv.y)
ht <- comb.fun(ht, "Hv")
summary(ht$Hv)  # 946 NAs

# Ks = Kl/Hv
names(ht)[names(ht) == "Ks"] <- "Ks.x"
ht$Ks.y <- (ht$Kl / ht$Hv) * 10000
summary(ht$Ks.x)
summary(ht$Ks.y)
ht <- comb.fun(ht, "Ks")
summary(ht$Ks)  #] 1164 NAs

# HSM completation
names(ht)[names(ht) == "HSM"] <- "HSM.x"
ht$HSM.y <- ht$MinWP_md - ht$P50
summary(ht$HSM.x)
summary(ht$HSM.y)
ht <- comb.fun(ht, "HSM")
summary(ht$HSM)  #] 1691 NAs

# Order
ht <- ht[, c("Species", "genus", "family", "order", "group", "Growth.form", "life.form", "leaf.form", "biome", "Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "AlAs", "WD", 
             "Vdia", "VD", "VLmax", "Hv", "HSM", "PItlp", "WUE", "SLA", "Aarea", "Amass", "LL", "pLL", "N", "Narea", "Rd", "gmax")]

write.csv(ht, "data/ht/HydraTRY_Pablo_traits.csv", row.names = F)
