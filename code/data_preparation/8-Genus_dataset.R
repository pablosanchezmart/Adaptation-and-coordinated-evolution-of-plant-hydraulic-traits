##################################### GENERA DATABASE #########################################################

### FUNCTIONS ####

## Mode function (for categorical variables)
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
# Maximum, minimum and mean function to avoid problems with NAs. 
max.fun <- function(x){
  y <- max(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}
min.fun <- function(x){
  y <- min(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}
mean.fun <- function(x){
  y <- mean(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}

#### Data ----------------------------------------------------------------------------------------------------- ####

ht <- read.csv("data/HydraTRY_Pablo_2020.csv", header = T)

### GENUS DATASET ####

# Aggregate data
catVar <- c("order", "family", "group", "Growth.form", "life.form", "leaf.form", "biome")
numVar <- names(ht)[c(10:length(ht))]

ht_aggr_mean <- aggregate(ht[, numVar], by = list(ht$genus), FUN = mean.fun)
ht_aggr_mode <- aggregate(ht[, catVar], by = list(ht$genus), FUN = mode.fun)

# Merge categorical and numerical variables
ht_aggr <- merge(ht_aggr_mean, ht_aggr_mode, by = "Group.1")

names(ht_aggr)[names(ht_aggr) == "Group.1"] <- "genus"
length(ht_aggr$genus) # 817

summary(ht_aggr)
# Order
ht_aggr <- ht_aggr[, names(ht)[-1]]
names(ht_aggr)

# Save
write.csv(ht_aggr, "data/HydraTRY_Genus_2020.csv", row.names = F)
print("====> data/HydraTRY_Genus_2020.csv")
