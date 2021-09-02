############################################ CLIMATE-HYDRAULIC FUNCTIONAL TRAITS PREVIOUS ANALYSES ###############################

# Pablo Sanchez Martínez
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

### MISSING DATA ####

# Total observations
totalObs <- length(ht[,1])

# Observations available for the whole potential predictors set

with(ht, length(which(!is.na(log_negP50) & is.na(PC1))))  # 83
with(ht, length(which(!is.na(log_negMinWP_md) & is.na(PC1))))  # 16
with(ht, length(which(!is.na(log_Ks) & is.na(PC1))))  # 50
with(ht, length(which(!is.na(log_Hv) & is.na(PC1))))  # 45

# Total observations
totalObs <- length(ht[, 1])
totalAng <- length(ht[which(ht$group == "Angiosperms"), 1])
totalGym <- length(ht[which(ht$group == "Gymnosperms"), 1])
totalPter <- length(ht[which(ht$group == "Pteridophytes"), 1])

total <- c(totalObs, totalAng, totalGym, totalPter)

# Observations available for the whole potential predictors set
mis.df <- data.frame()
n <- 0
for(i in Tvars){
  n <- n + 1
  completeObs_var <- length(which(!is.na(ht[, i])))
  Ang_obs <- length(which(!is.na(ht[ht$group == "Angiosperms", i])))
  Gym_obs <- length(which(!is.na(ht[ht$group == "Gymnosperms", i])))
  Pte_obs <- length(which(!is.na(ht[ht$group == "Pteridophytes", i])))
  mis.df <- rbind(mis.df, completeObs_var)
  mis.df[n, 2] <- Ang_obs
  mis.df[n, 3] <- Gym_obs
  mis.df[n, 4] <- Pte_obs
}
rownames(mis.df) <- Tvars
colnames(mis.df) <- c("Var Observations", "Angiosperms", "Gymnosperms", "Pterydophytes")
mis.df <- rbind(total, mis.df)
row.names(mis.df)[1] <- "Total Observations"
knitr::kable(mis.df)

### CORRELATION BETWEEN VARIABLES ####

### 1. Correlation between explanatory variables ###

# Create the correlation pairs plots

#HydVars
pairs(ht[, hydVars], lower.panel = corP, upper.panel = upperP)

#climVars
pairs(ht[, climVars], lower.panel = corP, upper.panel = upperP)

#soilVars
pairs(ht[, soilVars], lower.panel = corP, upper.panel = upperP)

#wholeset
pairs(ht[, Tvars], lower.panel = corP, upper.panel = upperP)

#analyses set
pairs(ht[, Avars], lower.panel = corP, upper.panel = upperP)

# Correlation matrix
corVars <- rcorr(as.matrix(ht[,Tvars]), type = "pearson")
# setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/")
# write.csv(corVars$r, "results1/Env-HT/var_correlation.csv")
# write.csv(corVars$n, "results1/Env-HT/var_correlation_n.csv")
# write.csv(corVars$p, "results1/Env-HT/var_correlation_p.csv")
# setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")


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
p

# gradient_color(c("#00AFBB", "#E7B800", "#FC4E07"))

## Correlations between principal components and hydraulic traits
hyd.coord[, 1:3]


### PRINCIPAL COMPONENTS REGRESSIONS ####

ht_pc <- read.csv("Ht_Pablo/HydraTRY_Pablo_2019.csv", header = T)

### log_Hmax ####

## Data
reg.Data <- ht_pc[, c("log_Hmax", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$log_Hmax, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "log_Hmax") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$log_Hmax, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "log_Hmax") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$log_Hmax, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "log_Hmax") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(log_Hmax ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$log_Hmax, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "log_Hmax", main = "log_Hmax", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$log_Hmax ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(log_Hmax ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "log_Hmax",
        ticktype = "detailed", color = "heat", main = "log_Hmax gam", se = 0)


### log_negP50 ####

## Data
reg.Data <- ht_pc[, c("log_negP50", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$log_negP50, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "log_negP50") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$log_negP50, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "log_negP50") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$log_negP50, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "log_negP50") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(log_negP50 ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$log_negP50, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "log_negP50", main = "log_negP50", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$log_negP50 ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(log_negP50 ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "log_negP50",
        ticktype = "detailed", color = "heat", main = "log_negP50 gam", se = 0)

### log_Ks ####

## Data
reg.Data <- ht_pc[, c("log_Ks", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$log_Ks, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "log_Ks") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$log_Ks, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "log_Ks") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$log_Ks, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "log_Ks") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(log_Ks ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$log_Ks, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "log_Ks", main = "log_Ks", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$log_Ks ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(log_Ks ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "log_Ks",
        ticktype = "detailed", color = "heat", main = "log_Ks gam", se = 0)


### log_Kl ####

## Data
reg.Data <- ht_pc[, c("log_Ks", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$log_Kl, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "log_Kl") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$log_Kl, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "log_Kl") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$log_Kl, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "log_Kl") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(log_Kl ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$log_Kl, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "log_Kl", main = "log_Kl", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$log_Kl ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(log_Kl ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "log_Kl",
        ticktype = "detailed", color = "heat", main = "log_Kl gam", se = 0)


### log_negPItlp ####

## Data
reg.Data <- ht_pc[, c("log_negPItlp", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$log_negPItlp, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "log_negPItlp") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$log_negPItlp, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "log_negPItlp") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$log_negPItlp, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "log_negPItlp") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(log_negPItlp ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$log_negPItlp, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "log_negPItlp", main = "log_negPItlp", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$log_negPItlp ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(log_negPItlp ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "log_negPItlp",
        ticktype = "detailed", color = "heat", main = "log_negPItlp gam", se = 0)

### log_negMinWP_md ####

## Data
reg.Data <- ht_pc[, c("log_negMinWP_md", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$log_negMinWP_md, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "log_negMinWP_md") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$log_negMinWP_md, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "log_negMinWP_md") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$log_negMinWP_md, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "log_negMinWP_md") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(log_negMinWP_md ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$log_negMinWP_md, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "log_negMinWP_md", main = "log_negMinWP_md", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$log_negMinWP_md ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(log_negMinWP_md ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "log_negMinWP_md",
        ticktype = "detailed", color = "heat", main = "log_negMinWP_md gam", se = 0)

### Rd ####

## Data
reg.Data <- ht_pc[, c("Rd", "group", "PC1", "PC2", "PC3")]

## Correlation with the principal components
# PC 1
pc1 <- ggplot(reg.Data, aes(x = ht_pc$PC1, y = ht_pc$Rd, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC1", y = "Rd") + theme_bw()

# PC 2
pc2 <- ggplot(reg.Data, aes(x = ht_pc$PC2, y = ht_pc$Rd, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC2", y = "Rd") + theme_bw()

# PC3
pc3 <- ggplot(reg.Data, aes(x = ht_pc$PC3, y = ht_pc$Rd, col = ht_pc$group)) + geom_point() + 
  geom_smooth(method="lm") + labs(x = "PC3", y = "Rd") + theme_bw()

# Plot
ggarrange(pc1, pc2, pc3)

## GLM
pcReg <- lm(Rd ~ (PC1 + PC2 + PC3) * group, data = ht_pc)
summary(pcReg)
plot(pcReg)

# 3D scatter plot with regression plane
colors <- c("#00AFBB",  "#FC4E07", "#07fc4e")
colors <- colors[as.numeric(ht_pc$group)]
s3d <- scatterplot3d(x = ht_pc$PC1, y = ht_pc$PC2, z = ht_pc$Rd, type = "h", color = colors, 
                     angle=55, pch = 16, xlab = "PC1", ylab = "PC2", zlab = "Rd", main = "Rd", 
                     highlight.3d = F, col.grid = "grey", lty.grid=par("lty"), box = T)
# Add regression plane
my.lm <- lm(ht_pc$Rd ~ ht_pc$PC1 + ht_pc$PC2)
s3d$plane3d(my.lm)

## gam model
gamMod <- gam(Rd ~ s(PC1, bs="ps") + s(PC2,bs="ps") + s(PC3,bs="ps"), method="REML",data=ht_pc, family=gaussian())
summary(gamMod)

plot(gamMod, shade=T, scale=0, all.terms=T)

# 3D surface response
vis.gam(gamMod, view = c("PC2", "PC1"), n.grid = 50, theta = 35, phi = 32, zlab = "Rd",
        ticktype = "detailed", color = "heat", main = "Rd gam", se = 0)
