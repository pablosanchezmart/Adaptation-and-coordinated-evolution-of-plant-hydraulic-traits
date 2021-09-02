#### INITIALIZATION ###################################################################################################################

# P Sanchez-Martinez

### PACKAGES -------------------------------------------------------------------------------------------------------####

library(plyr)
library(dplyr)
library(ape)
library(geiger)
library(MCMCglmm)
library(stringr)
library(tictoc)
library(ggplot2)
library(ggsci)
library(factoextra)
library(phytools)
library(doParallel)
library(foreach)
library(Hmisc)
library(scales)
library(ggsci)
library(plotbiomes)
library(ggpubr)
library(RColorBrewer)
library(V.PhyloMaker)
library(bayestestR)

print("Running initialization...")
options(scipen = 999, digits = 3)

#### FUNCTIONS ----------------------------------------------------------------------------------------------------- ####

# Complete data
completePhyloFun <- function(tree, dataset, desiredCols) {
  phyData <- list()
  completeVec <- complete.cases(dataset[, desiredCols])
  completeData <- dataset[completeVec, ]
  phylo <- ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% completeData$animal])
  phyData$dta <- as.data.frame(completeData[completeData$animal %in% phylo$tip.label, ])
  phyData$phylo <- phylo
  return(phyData)
}

#### DIRECTIORIES -------------------------------------------------------------------------------------------------- ####

# Species datset
ht <- read.csv("data/HydraTRY_Pablo_2020.csv", header = T)  # Includes three principal components

# Species level phylogeny 

mega.tr <- GBOTB.extended

# Genus random effect column
plmm.data <- ht
plmm.data$genus.Rand <- plmm.data$genus
plmm.data$animal <- plmm.data$genus
plmm.data$Spp <- str_replace_all(plmm.data$Species, " ", "_")

# Output directories
output.dir <- "outputs/"
manuscript.dir <- "manuscript/"

#### INTERSECTION DATA --------------------------------------------------------------------------------------------- ####

plmm.data <- plmm.data[plmm.data$Spp %in% mega.tr$tip.label, ]
ht.tree <- drop.tip(mega.tr, setdiff(mega.tr$tip.label, plmm.data$Spp))

## Prune at genus level

## get a list of all genera
tips <- ht.tree$tip.label
genera <- unique(sapply(strsplit(tips,"_"),function(x) x[1]))

## now drop all but one of each
ii<-sapply(genera,function(x,y) grep(x,y)[1],y=tips)
ht.tree <- drop.tip(ht.tree,setdiff(ht.tree$tip.label,tips[ii]))

ht.tree$tip.label<-sapply(strsplit(ht.tree$tip.label,"_"),function(x) x[1])

ht.tree$node.label <- NULL
length(plmm.data$Species) # 1356 rows
plmm.data$family <- NULL

#### MODEL SPECIFICATIONS ------------------------------------------------------------------------------------------- ####

### PRIORS

# Phylogenetic mixed multi models
priorMulti2 <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2), 
         G2 = list(V=diag(2)/2,nu=2)))

fmly <- "gaussian"
fmly2 <- c("gaussian", "gaussian")

# nitter = 8000000
# burni = 1000000
# thinni = 4000
# set.seed(123)

nitter = 4000000
burni = 100000
thinni = 2000
set.seed(123)

nitter = 200
burni = 100
thinni = 10
set.seed(123)

#### MODELS DESCRIPTION  --------------------------------------------------------------------------------------####

variables <- c("PC1", "PC2", "PC3", "log_Ks", "log_negP50", "log_Hv", "log_negMinWP_md", "HSM", "log_Kl")

# Bivariate models
multiMdls <- data.frame(rbind(
  cbind("Var_covariance", "NG", "log_Ks, log_negP50", 2, "cbind(log_Ks, log_negP50) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, log_Hv", 2, "cbind(log_Ks, log_Hv) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, log_negMinWP_md", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, PC1", 2, "cbind(log_Ks, PC1) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, PC1", 2, "cbind(log_Ks, PC1) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, PC2", 2, "cbind(log_Ks, PC2) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, PC2", 2, "cbind(log_Ks, PC2) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, PC3", 2, "cbind(log_Ks, PC3) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, PC3", 2, "cbind(log_Ks, PC3) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, log_Hv", 2, "cbind(log_negP50, log_Hv) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, log_negMinWP_md", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, PC1", 2, "cbind(log_negP50, PC1) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, PC1", 2, "cbind(log_negP50, PC1) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, PC2", 2, "cbind(log_negP50, PC2) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, PC2", 2, "cbind(log_negP50, PC2) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, PC3", 2, "cbind(log_negP50, PC3) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, PC3", 2, "cbind(log_negP50, PC3) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Hv, log_negMinWP_md", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Hv, PC1", 2, "cbind(log_Hv, PC1) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Hv, PC1", 2, "cbind(log_Hv, PC1) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Hv, PC2", 2, "cbind(log_Hv, PC2) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Hv, PC2", 2, "cbind(log_Hv, PC2) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Hv, PC3", 2, "cbind(log_Hv, PC3) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Hv, PC3", 2, "cbind(log_Hv, PC3) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negMinWP_md, PC1", 2, "cbind(log_negMinWP_md, PC1) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negMinWP_md, PC1", 2, "cbind(log_negMinWP_md, PC1) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negMinWP_md, PC2", 2, "cbind(log_negMinWP_md, PC2) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negMinWP_md, PC2", 2, "cbind(log_negMinWP_md, PC2) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negMinWP_md, PC3", 2, "cbind(log_negMinWP_md, PC3) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negMinWP_md, PC3", 2, "cbind(log_negMinWP_md, PC3) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Ks, log_negP50, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Ks, log_Hv, PC1, PC2, PC3", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Ks, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_negP50, log_Hv, PC1, PC2, PC3", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_negP50, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Hv, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand")
), stringsAsFactors =  T)
names(multiMdls) <- c("type", "group","variables", "N_ResVars", "fix.frml","ran.frml", "NP_ran.frml")
name <- paste0(as.character(multiMdls[, "type"]), "_", as.character(multiMdls[, "group"]), "_", as.character(multiMdls[, "variables"]))
name <- str_replace(name, ", ", "_")
multiMdls$name <- name
multiMdls

