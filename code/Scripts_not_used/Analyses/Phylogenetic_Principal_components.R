############################# PHYLOGENETIC PRINCIPAL COMPONENTS ANALYSES #############################

# Pablo Sanchez
# 05-2019

#### PACKAGES ####
library(phytools)
library(geiger)
#pca
library(factoextra)
library(ggbiplot)

load("RData/Phylogenetic_pca.RData")
### PARAMETERS #####

hydVars <- c("log_Hmax", "log_negP50", "log_Ks", "log_Kl", "log_Hv", "log_negPItlp", "log_negMinWP_md", "Rd", "log_gmax")

climVars <- c("sqrt_Prec.WarmQ" ,"TMax.warmest", "log_MAP", "log_Prec.DryQ", "MAT", "AI","srad", "windMax", "VPDMax", "TRange")

soilVars <- c("ADBRock", "ph.200cm", "clay.60cm_perc", "silt.60cm_perc",  "sand.60cm_perc", "SWC.200cm")

Tvars <- c("log_Hmax", "log_negP50", "log_Ks", "log_Kl", "log_Hv", "log_negPItlp", "log_negMinWP_md", "log_gmax", "Rd",                                                                        # Hydraulic traits
           "sqrt_Prec.WarmQ" ,"TMax.warmest", "log_MAP", "log_Prec.DryQ", "MAT", "AI","srad", "windMax", "VPDMax", "TRange",  # Clim vars
           "ADBRock", "ph.200cm", "clay.60cm_perc", "silt.60cm_perc",  "sand.60cm_perc", "SWC.200cm")                         # Soil vars

Avars <- c("group", "log_Hmax", "log_negP50", "log_Ks", "log_Hv", "log_negPItlp", "log_negMinWP_md", "log_gmax", "Rd",                         # Hyd vars                                                 # Hydraulic traits
           "sqrt_Prec.WarmQ" ,"TMax.warmest", "log_MAP", "log_Prec.DryQ", "MAT", "AI","srad", "windMax", "VPDMax", "TRange",  # Clim vars
           "ADBRock", "ph.200cm", "clay.60cm_perc", "silt.60cm_perc",  "sand.60cm_perc", "SWC.200cm")                         # Soil vars

envVars <- c("sqrt_Prec.WarmQ" ,"TMax.warmest", "log_MAP", "log_Prec.DryQ", "MAT", "AI","srad", "windMax", "VPDMax", "TRange",
             "ADBRock", "ph.200cm", "clay.60cm_perc", "silt.60cm_perc",  "sand.60cm_perc", "SWC.200cm")

### DATA ####
setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")

# Genus datset
ht <- read.csv("Ht_Pablo/HydraTRY_Genus.csv", header = T)  # Includes three principal components
length(ht$Genus)
names(ht)

## Data
na.obs <- which(is.na(ht$sqrt_Prec.WarmQ))
ht_pca <- ht[-na.obs, envVars]
ht_pca_hyd <- ht[-na.obs , hydVars]
ht_pca_gr <- ht$group[-na.obs]
ht_pca_de <- ht$leaf.form[-na.obs]
ht_pc <- ht[-na.obs, ]

na.obs <- which(is.na(ht_pca$AI))
ht_pca <- ht_pca[-na.obs, envVars]
ht_pca_hyd <- ht_pca_hyd[-na.obs , hydVars]
ht_pca_gr <- ht_pca_gr[-na.obs]
ht_pca_de <- ht_pca_de[-na.obs]
ht_pc <- ht_pc[-na.obs, ]

row.names(ht_pca) <- ht_pc$Genus

# Genus tree
ht.tree <- read.tree("Phylogeny_data/Ht_tree/hydratry_genus.tre") # Genus

# Intersenction data
ppca.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% row.names(ht_pca)])
ppca.data <- ht_pca[ppca.tree$tip.label, ]

# Check names of tree and ht
namesCheck <- name.check(ppca.tree, ppca.data)
# Taxa in tree but not in data
length(namesCheck$tree_not_data)
# Taxa in data but not in tree (NAs)
length(namesCheck$data_not_tree)

### PPCA ####
ppca.data <- as.matrix(ppca.data)

ppca <- phyl.pca(ppca.tree, ppca.data, method = "lambda", mode = "corr")

#Eigenvalues

summary(ppca)
plot(ppca)

# Variables correlation (loadings)
knitr::kable(ppca$L[, 1:3])

# Lambda
ppca$lambda
ppca$logL.lambda
# Plot
biplot(ppca)

# Genus scores
ppca$S

ppcaCoords <- as.data.frame(ppca$S[, 1:3])
ppcaCoords$genus <- rownames(ppcaCoords)

write.csv(ppcaCoords, "Outputs/ppcaScores.ssv", row.names = T)

# Join ppca scores with hydratry database

#save.image(file = "RData/Phylogenetic_pca.RData")
