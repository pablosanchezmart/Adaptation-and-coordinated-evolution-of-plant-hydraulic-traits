#### INITIALIZATION ###################################################################################################################

# P Sanchez-Martinez

### PACKAGES -------------------------------------------------------------------------------------------------------####

# devtools::install_github("valentinitnelav/plotbiomes")
# install.packages("bayestestR")

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
ht.gen <- read.csv("data/HydraTRY_Genus_2020.csv", header = T)

# Genus tree
ht.tree <- read.tree("data/hydratry_genus.tre") # Genus

# Output directories
output.dir <- "outputs/"
manuscript.dir <- "manuscript/"

#### DATABASE LENGTH AND COVERAGE ---------------------------------------------------------------------------------- ####

length(ht$Species)
length(unique(ht$genus))
length(unique(ht$family))

#### INTERSECTION DATA --------------------------------------------------------------------------------------------- ####

phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% ht$genus])
plmm.data <- ht[ht$genus %in% phylo$tip.label, ]
length(plmm.data$genus)
plmm.data$fmly <- plmm.data$family
plmm.data$family <- NULL # gives problems with MCMCglmm

# Genus random effect column
plmm.data$genus.Rand <- plmm.data$genus
plmm.data$animal <- plmm.data$genus

#### MODEL SPECIFICATIONS ------------------------------------------------------------------------------------------- ####

### PRIORS

# Non phylogenetic mixed models
prior1 <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 0.002)))

# Phylogenetic mixed models
prior2 <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 0.002), 
           G2 = list(V = 1, nu = 0.002))) # two random effects, phylogeny and genus

# Phylogenetic mixed multi models
priorMulti2 <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2), 
         G2 = list(V=diag(2)/2,nu=2)))



# Non-Phylogenetic mixed multi models
priorMulti2_NP <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2)))

## Ohter non informative priors

# With parameter expansion for random effects
# priorMulti2 <- list(
#   R=list(V=diag(2)/2,nu=2),
#   G=list(G1=list(V=diag(2)/2,nu=2, alpha.mu = rep(0,2), alpha.V = diag(25^2, 2, 2)),
#          G2 = list(V=diag(2)/2,nu=2, alpha.mu = rep(0,2), alpha.V = diag(25^2, 2, 2))))

# Cauchy
priorCau <- list(
  R = list(V = 1, nu = 0.002),
  G = list(G1 = list(V= 1, nu = 1, alpha.mu = 0),
           G2 = (G1 = list(V= 1, nu = 1, alpha.mu = 0))))


fmly <- "gaussian"
fmly2 <- c("gaussian", "gaussian")

nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

# nitter = 200
# burni = 100
# thinni = 10
# set.seed(123)


#### MODELS DESCRIPTION  --------------------------------------------------------------------------------------####

variables <- c("PC1", "PC2", "PC3", "log_Ks", "log_negP50", "log_Hv", "log_negMinWP_md", "HSM", "log_Kl")

# Univariate models 
uniMdls <- data.frame(rbind(
  cbind("Var_variance_patitioning", "NG", "PC1", 1, "PC1 ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "PC2", 1, "PC2 ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "PC3", 1, "PC3 ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "log_Ks", 1, "log_Ks ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "log_negP50", 1, "log_negP50 ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "log_Hv", 1, "log_Hv ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "log_negMinWP_md", 1,  "log_negMinWP_md ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "HSM", 1,  "HSM ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Var_variance_patitioning", "NG", "log_Kl", 1,  "log_Kl ~ 1", "~ animal + genus.Rand", "~ genus.Rand"),
  
  cbind("Ht_Env_1","NG", "log_Ks, PC1", 1, "log_Ks ~ PC1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Ks, PC2", 1, "log_Ks ~ PC2", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Ks, PC3", 1, "log_Ks ~ PC3", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_negP50, PC1", 1, "log_negP50 ~ PC1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_negP50, PC2", 1, "log_negP50 ~ PC2", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_negP50, PC3", 1, "log_negP50 ~ PC3", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Hv, PC1", 1, "log_Hv ~ PC1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Hv, PC2", 1, "log_Hv ~ PC2", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Hv, PC3", 1, "log_Hv ~ PC3", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Ks, PC1", 1, "log_Ks ~ PC1 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Ks, PC2", 1, "log_Ks ~ PC2 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Ks, PC3", 1, "log_Ks ~ PC3 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_negP50, PC1", 1, "log_negP50 ~ PC1 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_negP50, PC2", 1, "log_negP50 ~ PC2 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_negP50, PC3", 1, "log_negP50 ~ PC3 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Hv, PC1", 1, "log_Hv ~ PC1 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Hv, PC2", 1, "log_Hv ~ PC2 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Hv, PC3", 1, "log_Hv ~ PC3 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Ks, PC1", 1, "log_Ks ~ PC1 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Ks, PC2", 1, "log_Ks ~ PC2 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Ks, PC3", 1, "log_Ks ~ PC3 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_negP50, PC1", 1, "log_negP50 ~ PC1 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_negP50, PC2", 1, "log_negP50 ~ PC2 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_negP50, PC3", 1, "log_negP50 ~ PC3 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Hv, PC1", 1, "log_Hv ~ PC1 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Hv, PC2", 1, "log_Hv ~ PC2 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Hv, PC3", 1, "log_Hv ~ PC3 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "HSM, PC1", 1, "HSM ~ PC1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "HSM, PC2", 1, "HSM ~ PC2", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "HSM, PC3", 1, "HSM ~ PC3", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Kl, PC1", 1, "log_Kl ~ PC1", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Kl, PC2", 1, "log_Kl ~ PC2", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","NG", "log_Kl, PC3", 1, "log_Kl ~ PC3", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "HSM, PC1", 1, "HSM ~ PC1 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "HSM, PC2", 1, "HSM ~ PC2 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "HSM, PC3", 1, "HSM ~ PC3 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Kl, PC1", 1, "log_Kl ~ PC1 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Kl, PC2", 1, "log_Kl ~ PC2 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","G", "log_Kl, PC3", 1, "log_Kl ~ PC3 + group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "HSM, PC1", 1, "HSM ~ PC1 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "HSM, PC2", 1, "HSM ~ PC2 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "HSM, PC3", 1, "HSM ~ PC3 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Kl, PC1", 1, "log_Kl ~ PC1 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Kl, PC2", 1, "log_Kl ~ PC2 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Env_1","GI", "log_Kl, PC3", 1, "log_Kl ~ PC3 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  
  cbind("Ht_Ht_1","NG", "log_Ks, log_negP50", 1, "log_Ks ~ log_negP50", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Ht_1","GI", "log_Ks, log_negP50", 1, "log_Ks ~ log_negP50 * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Ht_1","NG", "log_Ks, log_Hv", 1, "log_Ks ~ log_Hv", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Ht_1","GI", "log_Ks, log_Hv", 1, "log_Ks ~ log_Hv * group", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Ht_1","NG", "log_negMinWP_md, log_negP50", 1, "log_negMinWP_md ~ log_negP50", "~ animal + genus.Rand", "~ genus.Rand"),
  cbind("Ht_Ht_1","GI", "log_negMinWP_md, log_negP50", 1, "log_negMinWP_md ~ log_negP50 * group", "~ animal + genus.Rand", "~ genus.Rand")
  ), stringsAsFactors =  T)
names(uniMdls) <- c("type", "group","variables", "N_ResVars", "fix.frml","ran.frml", "NP_ran.frml")
name <- paste0(as.character(uniMdls[, "type"]), "_", as.character(uniMdls[, "group"]), "_", as.character(uniMdls[, "variables"]))
name <- str_replace(name, ", ", "_")
uniMdls$name <- name

# Bivariate models
multiMdls <- data.frame(rbind(
  cbind("Var_covariance", "NG", "log_Ks, log_negP50", 2, "cbind(log_Ks, log_negP50) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, log_negP50", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, log_Hv", 2, "cbind(log_Ks, log_Hv) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, log_Hv", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, log_negMinWP_md", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, log_negMinWP_md", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, PC1", 2, "cbind(log_Ks, PC1) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, PC1", 2, "cbind(log_Ks, PC1) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, PC2", 2, "cbind(log_Ks, PC2) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, PC2", 2, "cbind(log_Ks, PC2) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Ks, PC3", 2, "cbind(log_Ks, PC3) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Ks, PC3", 2, "cbind(log_Ks, PC3) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, log_Hv", 2, "cbind(log_negP50, log_Hv) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, log_Hv", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, log_negMinWP_md", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, log_negMinWP_md", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, PC1", 2, "cbind(log_negP50, PC1) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, PC1", 2, "cbind(log_negP50, PC1) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, PC2", 2, "cbind(log_negP50, PC2) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, PC2", 2, "cbind(log_negP50, PC2) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_negP50, PC3", 2, "cbind(log_negP50, PC3) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_negP50, PC3", 2, "cbind(log_negP50, PC3) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "NG", "log_Hv, log_negMinWP_md", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance", "G", "log_Hv, log_negMinWP_md", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
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
  cbind("Var_covariance_env", "NG", "log_Ks, log_negP50, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "G", "log_Ks, log_negP50, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + trait:group + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Ks, log_negP50, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "NG", "log_Ks, log_Hv, PC1, PC2, PC3", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "G", "log_Ks, log_Hv, PC1, PC2, PC3", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + trait:group + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Ks, log_Hv, PC1, PC2, PC3", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "NG", "log_Ks, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "G", "log_Ks, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + trait:group + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Ks, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "NG", "log_negP50, log_Hv, PC1, PC2, PC3", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "G", "log_negP50, log_Hv, PC1, PC2, PC3", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + trait:group + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_negP50, log_Hv, PC1, PC2, PC3", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "NG", "log_negP50, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "G", "log_negP50, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + trait:group + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_negP50, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "NG", "log_Hv, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "G", "log_Hv, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + trait:group + trait:PC1 + trait:PC2 + trait:PC3", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand"),
  cbind("Var_covariance_env", "GI", "log_Hv, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal + us(trait):genus.Rand", "~ us(trait):genus.Rand")
), stringsAsFactors =  T)
names(multiMdls) <- c("type", "group","variables", "N_ResVars", "fix.frml","ran.frml", "NP_ran.frml")
name <- paste0(as.character(multiMdls[, "type"]), "_", as.character(multiMdls[, "group"]), "_", as.character(multiMdls[, "variables"]))
name <- str_replace(name, ", ", "_")
multiMdls$name <- name
multiMdls

