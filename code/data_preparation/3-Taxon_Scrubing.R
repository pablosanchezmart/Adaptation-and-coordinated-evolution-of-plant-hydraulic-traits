#### SPECIES NAMES CHECKING AND TAXON SCRUBBING ################################################################

# P Sanchez-Martinez

print("Taxon common nomenclature and scrubbing... ")

### PACKAGES ####
# install.packages("Taxonstand", dependencies = T)
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("wcornwell/taxonlookup")
# install.packages("taxize", dependencies = T)


library(Taxonstand) # step 1: clean up nomenclature
#https://cran.r-project.org/web/packages/Taxonstand/Taxonstand.pdf
#http://viktoriawagner.weebly.com/blog/cleaning-species-names-with-r-i-taxonstand
library(taxonlookup) # step 2: look it up
#https://github.com/traitecoevo/taxonlookup/blob/master/README.md
#library(phyndr) # step 3: swap species in tree for others with traits
#https://github.com/traitecoevo/phyndr-ms/blob/master/vignette/vignette.md

library(taxize) # OPTION 2: clean nomenclature
#http://viktoriawagner.weebly.com/blog/cleaning-species-names-with-r-ii-taxize

library(dplyr)
library(stringr)

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
mean.fun <- function(x){
  y <- mean(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}

# Combination function
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

### USING TAXONSTAND ####

ht <- read.csv("data/ht/HydraTRY_Pablo_traits.csv", header = T)  # Last version
length(ht$Species) # 2305

htSpp <- TPL(ht$Species) # Look for new names in The Plant List

# New names (corrected)
ht$newSpecies <- paste(htSpp$New.Genus, htSpp$New.Species, sep=" ") # Add new names to dataframe
dif <- ht %>% filter(ht$Species != ht$newSpecies) %>% select(Species, newSpecies)
dif
# Manual edit
ht[which(ht$Species == "Populus BAL"), "newSpecies"] <- "Populus balsamifera"
ht[which(ht$Species == "Populus FRI"), "newSpecies"] <- "Populus spp"
ht[which(ht$Species == "Populus GAV"), "newSpecies"] <- "Populus deltoides x nigra"
ht[which(ht$Species == "Populus GHO"), "newSpecies"] <- "Populus deltoides x nigra"
ht[which(ht$Species == "Salix BJO"), "newSpecies"] <- "Salix spp"
ht[which(ht$Species == "Salix DEL"), "newSpecies"] <- "Salix delnortensis"
ht[which(ht$Species == "Salix Q83"), "newSpecies"] <- "Salix triandra x viminalis"
ht[which(ht$Species == "Salix TOR"), "newSpecies"] <- "Salix babylonica"
ht[which(ht$Species == "Eleaegnus umbellata"), "newSpecies"] <- "Elaeagnus umbellata"
ht[which(ht$Species == "Oriza globerrima"), "newSpecies"] <- "Oryza globerrima"
ht[which(ht$Species == "Oriza sativa"), "newSpecies"] <- "Oryza sativa"

# Add new names
ht$Species <- ht$newSpecies
ht$genus <- htSpp$New.Genus                                      # Add new names to dataframe

ht_unique <- ht[!duplicated(ht$Species), ]
length(ht$Species) - length(ht_unique$Species) # 144 species duplicated

# Duplicates dataframe
duplSpp <- duplicated(ht$Species)
length(which(duplSpp == TRUE))
dupl.df <- ht[duplicated(ht$Species), ]
dupl.df <- arrange(dupl.df, Species)

# New families (avoiding empty new data)
ht$family <- as.character(ht$family)
for(i in 1:length(ht$family)){
  if(htSpp$Family[i] == ""){
    ht$family[i] <- ht$family[i]
  } else {
    ht$family[i] <- htSpp$Family[i]
  }
}
ht$family <- as.factor(ht$family)

# Aggregate data
names(ht)
catVar <- c("order", "family", "genus", "group", "Growth.form", "life.form", "leaf.form", "biome")
numVar <- names(ht)[c(10:34)]

ht_aggr_mean <- aggregate(ht[, numVar], by = list(ht$Species), FUN = mean.fun)
ht_aggr_mode <- aggregate(ht[, catVar], by = list(ht$Species), FUN = mode.fun)

# Merge categorical and numerical variables
ht_aggr <- merge(ht_aggr_mean, ht_aggr_mode, by = "Group.1")

names(ht_aggr)[names(ht_aggr) == "Group.1"] <- "Species"
length(ht_aggr$Species)

## Checking
head(ht_aggr)

ht[which(ht$Species == "Acer saccharum"), "MinWP_pd"]
ht_aggr[which(ht_aggr$Species == "Acer saccharum"),  "MinWP_pd"]

ht_tax <- arrange(ht_aggr, Species)
names(ht_tax)

### TAXONLOOKUP TO FILL TAXONOMIC GROUPS GAPS ---------------------------------------------------------------------- ####

taxVars <- c("genus", "family", "order", "group")

# Taxonlookup
species <- as.character(ht_tax$Species)
taxonomy <- lookup_table(species, by_species = TRUE)
taxonomy <- cbind("Species" = rownames(taxonomy), data.frame(taxonomy, row.names=NULL))
head(taxonomy)

# Character
for(i in taxVars){
  ht_tax[, i] <- as.character(ht_tax[, i])
}
# Merge
ht_tax_compl <- merge(ht_tax, taxonomy, by = "Species", all.x = T)

# Combinate
for(i in taxVars){
  ht_tax_compl <- comb.fun(ht_tax_compl, i)
}
# Factors
for(i in taxVars){
  ht_tax_compl[, i] <- as.factor(ht_tax_compl[, i])
}

# Checking
summary(ht_tax_compl$genus) # 0 NAs
summary(ht_tax_compl$family) # 0 NAs
summary(ht_tax_compl$order) # 7 NAs
summary(ht_tax_compl$group) # 0 NAs
head(ht_tax_compl$genus)

# Order
ht_tax_compl <- ht_tax_compl[, c("Species", "genus", "family", "order", "group", "Growth.form", "life.form", "leaf.form", "biome", "Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "AlAs", "WD", 
                     "Vdia", "VD", "VLmax", "Hv", "HSM", "PItlp", "WUE", "SLA", "Aarea", "Amass", "LL", "pLL", "N", "Narea", "Rd", "gmax")]


#### EXCLUDE SOME SPECIES ---------------------------------------------------------------------------------------------------- ####

ht <- ht_tax_compl

# Exclude pteridophytes
ht$group <- as.character(ht$group)
ht_pter <- ht[which(ht$group != "Pteridophytes"), ]
ht_pter$group <- as.factor(ht_pter$group)
summary(ht_pter$group)

# Exclude non-woody
ht_pter$Growth.form <- as.factor(ht_pter$Growth.form)
ht_pter$life.form <- as.factor(ht_pter$life.form)

summary(ht_pter$Growth.form)
summary(ht_pter$life.form)

# Consider all trees and shrubs as woody
ht_pter$Growth.form[which(!is.na(ht_pter$life.form))] <- "W"
summary(ht_pter$Growth.form)

ht_pter_w <- ht_pter[which(ht_pter$Growth.form == "W"), ]
summary(ht_pter_w$Growth.form)
head(ht_pter_w)

# Save
write.csv(ht_pter_w, "data/ht/HydraTRY_Pablo_traits_taxscrub.csv", row.names = F)



