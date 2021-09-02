#### PRINCIPAL COMPONENTS CALCULATION ###################################################################################################################

# P Sanchez-Martinez

print("Running PCA computation...")

### FUNCTIONS -----------------------------------------------------------------------------------------------------------#### 

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}


### VARIABLES AND DATA --------------------------------------------------------------------------------------------------####

ht <- read.csv("data/HydraTRY_Pablo_2020.csv", header = T)

hydVars <- c("log_negP50", "log_Ks", "log_Kl", "log_Hv", "log_negMinWP_md", "HSM")   
envVars <- c("MAT","TMean_Warmest","TMean_Coldest", "AP", "log_Prec_Wettest","sqrt_Prec_Driest", "log_Prec_Seasonality", 
             "log_Prec_WetQ", "sqrt_Prec_DryQ", "sqrt_Prec_WarmQ", 
             "sqrt_Prec_ColdQ", "Diurnal_TRange", "Isotermality", "log_TSeasonality", "TMax.Warmest", "TMin.Coldest", "log_TRange",
             "TMean_Wettest", "TMean_Driest", "srad", "sradMin", "tmax", "tmin", "vapr", "windMax", "wind", "VPDMax", 
             "ABDRock", "clay_60cm", "ph_60cm", "sand_60cm", "SWC_200cm", "AI") 

envVars <- c("sqrt_Prec_WarmQ","log_Prec_Wettest", "tmax", "log_TSeasonality", "AP", "sqrt_Prec_DryQ", 
             "MAT", "AI","srad", "windMax", "VPDMax",
             "ABDRock", "ph_60cm", "clay_60cm", "sand_60cm", "SWC_200cm")

## Complete data
ht_pca <- completeFun(ht, envVars)
length(ht_pca$Species)

# Auxiliar database used later to merge PCs with HydraTRY
ht_pc <- ht_pca

# Environmental variables
ht_pca <- ht_pc[, envVars]

# Hydraulic traits
ht_pca_hyd <- ht_pc[, hydVars]

# For group plot
ht_pca_gr <- ht_pc$group
# For leaf form plot
ht_pca_de <- ht_pc$leaf.form

### PCA results ----------------------------------------------------------------------------------------------------####

pca <- prcomp(ht_pca, scale = T)
summary(pca)


#Eighenvalues
get_eig(pca)
fviz_eig(pca, addlabels = T, hjust=0.3)
dev.off()

# Results for Variables
res.var <- get_pca_var(pca)
round(res.var$cor[, 1:3] *-1, 3)            # Coordinates
round(res.var$contrib[, 1:3], 3)        # Contributions to the PCs
round(res.var$cos2, 3)           # Quality of representation 

# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
res.ind$cos2           # Quality of representation


### HT - PCs PEARSON CORRELATIONS ------------------------------------------------------------------------------------------------ ####

## Predict coordinates and compute cos2
hyd.coord <- cor(ht_pca_hyd, pca$x, use = "pairwise.complete.obs")
hyd.cos2 <- hyd.coord^2

## Correlations
hyd.coord[, 1:3]


### MERGE WITH HydraTRY DATASET ------------------------------------------------------------------------------------------#### 

## Save 3 PC

length(pca$x[, 1])
length(ht_pc$Species)
ht_pc$PC1 <- pca$x[, 1] * -1
ht_pc$PC2 <- pca$x[, 2] * -1
ht_pc$PC3 <- pca$x[, 3] * -1
summary(ht_pc)
length(ht_pc$Species)

## Join with Ht dataset

ht_pca <- ht_pc[, c("Species","PC1", "PC2", "PC3")]
ht_pca <- join(ht, ht_pca, by = "Species")
summary(ht_pca)
length(ht_pca$Species)

write.csv(ht_pca, "data/HydraTRY_Pablo_2020.csv", row.names = F)
print("====> data/HydraTRY_Pablo_2019.csv")

#### SAVE IMAGE ---------------------------------------------------------------------------------------------------------####

save.image(file = "data/RData/PCA.RData")
print("====> data/RData/PCA.RData")
