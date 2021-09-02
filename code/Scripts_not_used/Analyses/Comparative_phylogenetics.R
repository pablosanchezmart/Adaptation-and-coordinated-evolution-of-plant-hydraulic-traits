################################# PHYLOGENY ANALYSES ###################################################################################################


# Pablo Sanchez


### PACKAGES ####
#install.packages("caper", dependencies = T)
# devtools::install_github("richfitz/phyndr")
# devtools::install_github("ropenscilabs/datastorr")
# devtools::install_github("wcornwell/taxonlookup")
# install.packages("BiocManager")
# BiocManager::install("tidytree")
# BiocManager::install("ggtree")
# install.packages("phylosignal")

# tree
library(Taxonstand) # step 1: clean up nomenclature
library(taxonlookup) # step 2: look it up
library(phyndr) # step 3: swap species in tree for others with traits
library(taxize) # OPTION 2: clean nomenclature

library(dplyr)
#pgls
library(ape) # Phylogenetical analyses
library(geiger)
library(nlme)
library(phytools) # phylogenetical tools (force.ultrametric)
library(caper)

#phylogenetic signal
library(phylosignal)
library(adephylo)
library(phylobase)

# plotting
library(tidytree)
library(ggtree)
library(ggplot2)
library(ggpubr)


### FUNCTIONS ####

#Function to get complete cases:return subset without NA
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
# ggtree plot
treePlotFun <- function(tree, taxGroup, title){
  taxInfo <- lookup_table(tree$tip.label, by_species = TRUE)
  taxInfo[, taxGroup] <- as.factor(taxInfo[, taxGroup]) 
  groupInfo <- split(tree$tip.label, taxInfo[, taxGroup])
  tree.plot <- groupOTU(tree, groupInfo)
  tp <- ggtree(tree.plot, aes(colour = group), layout = "fan")
  tp <- tp + ggtitle(paste(title, sep = " ")) + theme(legend.position = "none") 
  plot(tp)
}

# Phylogenetic signal
phySigFun <- function(dta, txn = "Species", vrble, mthds){
  # Data
  ht.ps <- dta[, c(txn, paste(vrble))]
  ht.ps <- completeFun(ht.ps, c(txn, paste(vrble)))
  # Tree
  ht.psTree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% ht.ps[, txn]])
  # Order dataframe in the same pruned tree order
  ht.ps <- ht.ps[ht.psTree$tip.label, ]
  ht.ps <- ht.ps[, vrble]
  n <- length(ht.ps)
  # phylogntic signal
  phy.ps <- phylo4d(ht.psTree, ht.ps)
  phySig <- phyloSignal(phy.ps, methods = mthds)
  phySig$n <- n
  phySig <- as.data.frame(phySig)
  phySig <- round(phySig, digits = 3)
  phySig <- cbind(variable, phySig)
  return(phySig)
}


### PARAMETERS ####
hydVars <- c("log_Hmax", "log_negP50", "log_Ks","log_Kl", "log_Hv", "log_negPItlp", "log_negMinWP_md", "Rd", "log_gmax")
tVars <- c("log_Hmax", "log_negP50", "log_Ks","log_Kl", "log_Hv", "log_negPItlp", "log_negMinWP_md", "Rd", "log_gmax", "PC1", "PC2", "PC3")
### DATA ####

setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")
# Genus dataset
ht <- read.csv("Ht_Pablo/HydraTRY_Genus.csv", header = T)

row.names(ht) <- ht$Genus
ht.tree <- read.tree("Phylogeny_data/Ht_tree/hydratry_genus.tre") # Genus

treePlotFun(ht.tree, "family", "Phylogeny for hydraulics")

###PHYLOGENETIC SIGNAL OF INDIVIDUAL TRAITS ####

vars.ps <- data.frame()
for(variable in tVars){
 var.ps <- phySigFun(ht, txn = "Genus", variable, c("K", "Lambda"))
 vars.ps <- bind_rows(vars.ps, var.ps)
}
vars.ps

### CLIMATE PHYLIGENETICAL ANALYSES ####

# Check names of tree and ht
namesCheck <- name.check(ht.tree, ht)
# Taxa in tree but not in data
length(namesCheck$tree_not_data)
# Taxa in data but not in tree (NAs)
length(namesCheck$data_not_tree)


### log_Ks ####

## Data
gls.data <- completeFun(ht, c("log_Ks", "PC1"))
length(gls.data$log_Ks)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_Ks ~ PC1 + PC2 + PC3, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.5, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_Ks ~ PC1 + PC2 + PC3, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_Ks ~ PC1 + PC2 + PC3, correlation = pg.ht, data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = PC1, y = log_Ks, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC1, colour = model), data=coefs)
p2 <- scatterplot <- qplot(x = PC2, y = log_Ks, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC2, colour = model), data=coefs)
p3 <- scatterplot <- qplot(x = PC3, y = log_Ks, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC3, colour = model), data=coefs)

ggarrange(p1, p2, p3)


### log_Hv ####

## Data
gls.data <- completeFun(ht, c("log_Hv", "PC1"))
length(gls.data$log_Hv)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_Hv ~ PC1 + PC2 + PC3, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.9, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_Hv ~ PC1 + PC2 + PC3, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_Hv ~ PC1 + PC2 + PC3, correlation = pg.ht, 
              data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = PC1, y = log_Hv, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC1, colour = model), data=coefs)
p2 <- scatterplot <- qplot(x = PC2, y = log_Hv, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC2, colour = model), data=coefs)
p3 <- scatterplot <- qplot(x = PC3, y = log_Hv, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC3, colour = model), data=coefs)

ggarrange(p1, p2, p3)


### log_negP50 ####

## Data
gls.data <- completeFun(ht, c("log_negP50", "PC1"))
length(gls.data$log_negP50)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_negP50 ~ PC1 + PC2 + PC3, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.5, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_negP50 ~ PC1 + PC2 + PC3, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_negP50 ~ PC1 + PC2 + PC3, correlation = pg.ht, data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = PC1, y = log_negP50, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC1, colour = model), data=coefs)
p2 <- scatterplot <- qplot(x = PC2, y = log_negP50, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC2, colour = model), data=coefs)
p3 <- scatterplot <- qplot(x = PC3, y = log_negP50, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC3, colour = model), data=coefs)

ggarrange(p1, p2, p3)


### log_negMinWP_md ####

## Data
gls.data <- completeFun(ht, c("log_negMinWP_md", "PC1"))
length(gls.data$log_negMinWP_md)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_negMinWP_md ~ PC1 + PC2 + PC3, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.5, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_negMinWP_md ~ PC1 + PC2 + PC3, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_negMinWP_md ~ PC1 + PC2 + PC3, correlation = pg.ht, data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = PC1, y = log_negMinWP_md, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC1, colour = model), data=coefs)
p2 <- scatterplot <- qplot(x = PC2, y = log_negMinWP_md, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC2, colour = model), data=coefs)
p3 <- scatterplot <- qplot(x = PC3, y = log_negMinWP_md, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = PC3, colour = model), data=coefs)

ggarrange(p1, p2, p3)


### PGLS BETWEEN HYDRAULIC TRAITS ####

### P50 ~ Ks

## Data
gls.data <- completeFun(ht, c("log_negP50", "log_Ks"))
length(gls.data$log_Ks)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_Ks ~ log_negP50, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.5, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_Ks ~ log_negP50, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_Ks ~ log_negP50, correlation = pg.ht, data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = log_negP50, y = log_Ks, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = log_negP50, colour = model), data=coefs)
p1

### Hv ~ Ks

## Data
gls.data <- completeFun(ht, c("log_Hv", "log_Ks"))
length(gls.data$log_Ks)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_Hv ~ log_Ks, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.5, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_Hv ~ log_Ks, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_Hv ~ log_Ks, correlation = pg.ht, data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = log_Ks, y = log_Hv, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = log_Ks, colour = model), data=coefs)
p1

### MinWP ~ P50

## Data
gls.data <- completeFun(ht, c("log_negMinWP_md", "log_negP50"))
length(gls.data$log_negMinWP_md)
# Tree
gls.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% gls.data$Genus])
length(gls.tree$tip.label) 
gls.data <- gls.data[gls.tree$tip.label, ]
length(gls.data$Genus)

# Tree plot
#treePlotFun(gls.tree, "group", "Phylogeny for hydraulics")

### GLS NLME

#ols
ols <- lm(log_negP50 ~ log_negMinWP_md, data = gls.data)
summary(ols)

## covarince matrix gls
# Brownian motion (convergence)
bm.ht <- corBrownian(phy = gls.tree)
# Pagel lambda (intermediate)
pg.ht <- corPagel(0.5, phy = gls.tree, fixed = FALSE) # fixed = FALSE glm estimates lambda
#ou.ht <- corMartins(1, phy = gls.tree)

# Brwnian motion
bm.gls <- gls(log_negP50 ~ log_negMinWP_md, correlation = bm.ht, data = gls.data, method = "ML")
summary(bm.gls)

# Pagel independence
pg.gls <- gls(log_negP50 ~ log_negMinWP_md, correlation = pg.ht, data = gls.data, method = "ML")
summary(pg.gls)

# log-likelihood ratio tests
anova(bm.gls, ols)
anova(pg.gls, ols)
anova(bm.gls, pg.gls)

## Plot differences

## Coefficients
# bmPGLS
bmPgls.Coef <- data.frame(t(coef(bm.gls)))
# colnames(bmPgls.Coef) <- c("intercept", "slopePC1")
row.names(bmPgls.Coef) <- "bmPGLS"
# pgPGLS
pgPgls.Coef <- data.frame(t(coef(pg.gls)))
# colnames(pgPgls.Coef) <- c("intercept", "slopePC1")
row.names(pgPgls.Coef) <- "pgPGLS"
# OLS
ols.Coef <- data.frame(t(coef(ols)))
# colnames(ols.Coef) <- c("intercept", "slope")
row.names(ols.Coef) <- "OLS"

coefs <- rbind(bmPgls.Coef, pgPgls.Coef, ols.Coef)
coefs$model <- row.names(coefs)
coefs

## Plot
p1 <- scatterplot <- qplot(x = log_negMinWP_md, y = log_negP50, data = gls.data) + 
  geom_abline(aes(intercept = X.Intercept., slope = log_negMinWP_md, colour = model), data=coefs)
p1


### PHYLOGENETIC PATH ANALYSIS####


