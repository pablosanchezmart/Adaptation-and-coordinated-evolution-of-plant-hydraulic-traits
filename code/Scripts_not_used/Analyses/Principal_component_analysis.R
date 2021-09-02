############################################ PRINCIPAL COMPOENTENTS ANALYSIS ###############################

# Pablo Sanchez Mart?nez
# 14/03/2019

### PACKAGES ####
library(ggplot2)
library(nortest)
library(dplyr)
library(ggpubr)
library(mgcv)
library(Hmisc) # rcorr
#pca
library(factoextra)
# Could be interesting to plot pca
library(devtools)
# install_github("vqv/ggbiplot")
# join
library(plyr)
# Surface plot
library(scatterplot3d)
library(lattice)
library(plotly)

load("data/RData/PCA.RData")

### FUNCTIONS ####

# Correlation panel
corP <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, use = "pairwise.complete.obs"), digits=2)
  txt <- paste0("R = ", r)
  if(r > 0.5 | r < -0.5){
    text(0.5, 0.5, txt, col = "Red")
  } else {
    text(0.5, 0.5, txt)
  }
}
# Scatter plot panel
upperP <- function(x, y){
  points(x, y, pch = 19)
  # ggplot(ht, aes(x, y)) + geom_point + geom_smooth(method = "lm") + theme_bw()
}


### DATA ####
setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")

ht <- read.csv("Ht_Pablo/HydraTRY_Pablo_2019.csv", header = T)

# Variables to be used

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


### PRINCIPAL COMPONENT ANALYSES: CLIMATE+SOIL (environmental variables) ####

#http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

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

## PCA computation
pca <- prcomp(ht_pca, scale = T)
summary(pca)

# Save 3 PC
length(pca$x[, 1])
length(ht_pc$Species)
ht_pc$PC1 <- pca$x[, 1]
ht_pc$PC2 <- pca$x[, 2]
ht_pc$PC3 <- pca$x[, 3]
summary(ht_pc)
length(ht_pc$Species)

# Join with Ht dataset
ht_pca <- ht_pc[, c("Species","PC1", "PC2", "PC3")]
ht_pca <- join(ht, ht_pca, by = "Species")
summary(ht_pca)
length(ht_pca$Species)

#write.csv(ht_pca, "Ht_Pablo/HydraTRY_Pablo_2019.csv", row.names = F)

## PCA results
#Eighenvalues
get_eig(pca)
fviz_eig(pca, addlabels = T, hjust=0.3)
#dev.off()
# Results for Variables
res.var <- get_pca_var(pca)
res.var$cor            # Coordinates
res.var$contrib        # Contributions to the PCs
res.var$cos2           # Quality of representation 

# Save results
# setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Outputs")
# write.csv(res.var$cor, "PCA_variables_correlations.csv", row.names = F)
# write.csv(res.var$contrib, "PCA_variables_contributions.csv", row.names = F)

# Results for individuals
res.ind <- get_pca_ind(pca)
res.ind$coord          # Coordinates
res.ind$contrib        # Contributions to the PCs
# res.ind$cos2           # Quality of representation 

# PCA variables plot
fviz_pca_var(pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# ggbiplot(pca, obs.scale = 1, var.scale = 1,
#          groups = ht_pca_gr, ellipse = TRUE, circle = TRUE) +
#   scale_color_discrete(name = '') +
#   theme(legend.direction = 'horizontal', legend.position = 'top')

# biplot of individuals and variables
fviz_pca_biplot(pca,
                col.var = "#000000", col.ind = ht_pca_gr, geom = "point", palette = c("#00AFBB",  "#FC4E07", "#07fc4e"),
                addEllipses = TRUE, ellipse.level=0.95, 
                legend.title = "Groups",
                repel = T
)
# biplot including leaf form 
fviz_pca_biplot(pca,
                col.var = "#000000", col.ind = ht_pca_de, geom = "point", palette = c("#e74c3c", "#58c890", "#594ba0"),
                addEllipses = TRUE, ellipse.level=0.95, 
                legend.title = "Groups", 
                repel = T
)

### ADDING SUPPLEMENTARY VARIABLES TO SHOW CORRELATION WITH THE ALREADY COMPUTED PRINCIPAL COMPONENTS ####

## Predict coordinates and compute cos2
hyd.coord <- cor(ht_pca_hyd, pca$x, use = "pairwise.complete.obs")
hyd.cos2 <- hyd.coord^2

## Graph of variables including supplementary variables
p <- fviz_pca_ind(pca, geom = "point", col.ind = ht_pca_de, palette = c("#e74c3c", "#58c890", "#594ba0")
)
p <- fviz_add(p, hyd.coord*10, geom="arrow", repel = T, linetype = "solid", color = "black"
)

# gradient_color(c("#00AFBB", "#E7B800", "#FC4E07"))

## Correlations between principal components and hydraulic traits
hyd.coord[, 1:3]

#### SAVE IMAGE ####

save.image(file = "RData/PCA.RData")
