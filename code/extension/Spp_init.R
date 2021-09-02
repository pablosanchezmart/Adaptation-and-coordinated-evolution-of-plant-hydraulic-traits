#### INITIALIZATION ###################################################################################################################

# P Sanchez-Martinez

### PACKAGES -------------------------------------------------------------------------------------------------------####

# devtools::install_github("jinyizju/V.PhyloMaker")
library(plyr)
library(dplyr)
library(ape)
library(geiger)
library(MCMCglmm)
library(stringr)
library(tictoc)
library(V.PhyloMaker)

print("Running initialization...")

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

plmm.data <- ht
plmm.data$animal <- str_replace_all(plmm.data$Species, " ", "_")

# Output directories
output.dir <- "outputs/"
manuscript.dir <- "manuscript/"

#### INTERSECTION DATA --------------------------------------------------------------------------------------------- ####

plmm.data <- plmm.data[plmm.data$animal %in% mega.tr$tip.label, ]
ht.tree <- drop.tip(mega.tr, setdiff(mega.tr$tip.label, plmm.data$animal))
ht.tree$node.label <- NULL
length(plmm.data$Species) # 1356 rows
plmm.data$family <- NULL

### Join Agiosperms/gymnosperms and leaf form ###

plmm.data$group_leaf <- paste0(plmm.data$group, "_", plmm.data$leaf.form)
for(i in 1:length(plmm.data$group_leaf)){
  if(!is.na(str_extract(plmm.data$group_leaf[i], "NA"))){
    plmm.data$group_leaf[i] <- NA
  }
}
plmm.data$group_leaf <- as.character(plmm.data$group_leaf)
plmm.data[which(plmm.data$group_leaf == "Gymnosperms_D"), "group_leaf"] <- "Gymnosperms"
plmm.data[which(plmm.data$group_leaf == "Gymnosperms_E"), "group_leaf"] <- "Gymnosperms"
plmm.data[which(plmm.data$group_leaf == "Gymnosperms_W"), "group_leaf"] <- "Gymnosperms"
plmm.data$group_leaf <- as.factor(plmm.data$group_leaf)
summary(plmm.data$group_leaf)

#### MODEL SPECIFICATIONS ------------------------------------------------------------------------------------------- ####

### PRIORS

# Phylogenetic mixed models
prior2 <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 0.002))) # one random effects, phylogeny

# Non-Phylogenetic mixed models
prior2_NP <- list(
  R = list(V = 1, nu = 0.002)
) # one random effects, phylogeny

# Phylogenetic mixed multi models
priorMulti2 <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2)))

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

# nitter = 200
# burni = 100
# thinni = 10
# set.seed(123)

#### MODELS DESCRIPTION  --------------------------------------------------------------------------------------####

variables <- c("PC1", "PC2", "PC3", "log_Ks", "log_negP50", "log_Hv", "log_negMinWP_md", "HSM", "log_Kl")

# Univariate models 
Sp_uniMdls <- data.frame(rbind(
  cbind("Var_variance_patitioning", "NG", "PC1", 1, "PC1 ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "PC2", 1, "PC2 ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "PC3", 1, "PC3 ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "log_Ks", 1, "log_Ks ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "log_negP50", 1, "log_negP50 ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "log_Hv", 1, "log_Hv ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "log_negMinWP_md", 1,  "log_negMinWP_md ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "HSM", 1,  "HSM ~ 1", "~ animal", "_"),
  cbind("Var_variance_patitioning", "NG", "log_Kl", 1,  "log_Kl ~ 1", "~ animal", "_"),

  cbind("Ht_Env_1","NG", "log_Ks, PC1", 1, "log_Ks ~ PC1", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Ks, PC2", 1, "log_Ks ~ PC2", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Ks, PC3", 1, "log_Ks ~ PC3", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_negP50, PC1", 1, "log_negP50 ~ PC1", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_negP50, PC2", 1, "log_negP50 ~ PC2", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_negP50, PC3", 1, "log_negP50 ~ PC3", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Hv, PC1", 1, "log_Hv ~ PC1", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Hv, PC2", 1, "log_Hv ~ PC2", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Hv, PC3", 1, "log_Hv ~ PC3", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "HSM, PC1", 1, "HSM ~ PC1", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "HSM, PC2", 1, "HSM ~ PC2", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "HSM, PC3", 1, "HSM ~ PC3", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Kl, PC1", 1, "log_Kl ~ PC1", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Kl, PC2", 1, "log_Kl ~ PC2", "~ animal", "_"),
  cbind("Ht_Env_1","NG", "log_Kl, PC3", 1, "log_Kl ~ PC3", "~ animal", "_"),

  cbind("Ht_Env_1","G", "log_Ks, PC1", 1, "log_Ks ~ PC1 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Ks, PC2", 1, "log_Ks ~ PC2 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Ks, PC3", 1, "log_Ks ~ PC3 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_negP50, PC1", 1, "log_negP50 ~ PC1 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_negP50, PC2", 1, "log_negP50 ~ PC2 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_negP50, PC3", 1, "log_negP50 ~ PC3 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Hv, PC1", 1, "log_Hv ~ PC1 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Hv, PC2", 1, "log_Hv ~ PC2 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Hv, PC3", 1, "log_Hv ~ PC3 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "HSM, PC1", 1, "HSM ~ PC1 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "HSM, PC2", 1, "HSM ~ PC2 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "HSM, PC3", 1, "HSM ~ PC3 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Kl, PC1", 1, "log_Kl ~ PC1 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Kl, PC2", 1, "log_Kl ~ PC2 + group", "~ animal", "_"),
  cbind("Ht_Env_1","G", "log_Kl, PC3", 1, "log_Kl ~ PC3 + group", "~ animal", "_"),

  cbind("Ht_Env_1","GI", "log_Ks, PC1", 1, "log_Ks ~ PC1 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Ks, PC2", 1, "log_Ks ~ PC2 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Ks, PC3", 1, "log_Ks ~ PC3 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_negP50, PC1", 1, "log_negP50 ~ PC1 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_negP50, PC2", 1, "log_negP50 ~ PC2 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_negP50, PC3", 1, "log_negP50 ~ PC3 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Hv, PC1", 1, "log_Hv ~ PC1 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Hv, PC2", 1, "log_Hv ~ PC2 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Hv, PC3", 1, "log_Hv ~ PC3 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "HSM, PC1", 1, "HSM ~ PC1 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "HSM, PC2", 1, "HSM ~ PC2 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "HSM, PC3", 1, "HSM ~ PC3 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Kl, PC1", 1, "log_Kl ~ PC1 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Kl, PC2", 1, "log_Kl ~ PC2 * group", "~ animal", "_"),
  cbind("Ht_Env_1","GI", "log_Kl, PC3", 1, "log_Kl ~ PC3 * group", "~ animal", "_"),
  
  cbind("Ht_Env_1","GLI", "log_Ks, PC1", 1, "log_Ks ~ PC1 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Ks, PC2", 1, "log_Ks ~ PC2 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Ks, PC3", 1, "log_Ks ~ PC3 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_negP50, PC1", 1, "log_negP50 ~ PC1 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_negP50, PC2", 1, "log_negP50 ~ PC2 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_negP50, PC3", 1, "log_negP50 ~ PC3 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Hv, PC1", 1, "log_Hv ~ PC1 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Hv, PC2", 1, "log_Hv ~ PC2 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Hv, PC3", 1, "log_Hv ~ PC3 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_negMinWP_md, PC1", 1, "log_negMinWP_md ~ PC1 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_negMinWP_md, PC2", 1, "log_negMinWP_md ~ PC2 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_negMinWP_md, PC3", 1, "log_negMinWP_md ~ PC3 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "HSM, PC1", 1, "HSM ~ PC1 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "HSM, PC2", 1, "HSM ~ PC2 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "HSM, PC3", 1, "HSM ~ PC3 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Kl, PC1", 1, "log_Kl ~ PC1 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Kl, PC2", 1, "log_Kl ~ PC2 * group_leaf", "~ animal", "_"),
  cbind("Ht_Env_1","GLI", "log_Kl, PC3", 1, "log_Kl ~ PC3 * group_leaf", "~ animal", "_")
), stringsAsFactors =  T)
names(Sp_uniMdls) <- c("type", "group","variables", "N_ResVars", "fix.frml","ran.frml", "NP_ran.frml")
name <- paste0(as.character(Sp_uniMdls[, "type"]), "_", as.character(Sp_uniMdls[, "group"]), "_", as.character(Sp_uniMdls[, "variables"]))
name <- str_replace(name, ", ", "_")
Sp_uniMdls$name <- name

# Bivariate models
Sp_multiMdls <- data.frame(rbind(
  cbind("Sp_Var_covariance", "NG", "log_Ks, log_negP50", 2, "cbind(log_Ks, log_negP50) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Ks, log_Hv", 2, "cbind(log_Ks, log_Hv) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Ks, log_negMinWP_md", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Ks, PC1", 2, "cbind(log_Ks, PC1) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Ks, PC2", 2, "cbind(log_Ks, PC2) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Ks, PC3", 2, "cbind(log_Ks, PC3) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negP50, log_Hv", 2, "cbind(log_negP50, log_Hv) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negP50, log_negMinWP_md", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negP50, PC1", 2, "cbind(log_negP50, PC1) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negP50, PC2", 2, "cbind(log_negP50, PC2) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negP50, PC3", 2, "cbind(log_negP50, PC3) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Hv, log_negMinWP_md", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Hv, PC1", 2, "cbind(log_Hv, PC1) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Hv, PC2", 2, "cbind(log_Hv, PC2) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_Hv, PC3", 2, "cbind(log_Hv, PC3) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negMinWP_md, PC1", 2, "cbind(log_negMinWP_md, PC1) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negMinWP_md, PC2", 2, "cbind(log_negMinWP_md, PC2) ~ trait-1", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "NG", "log_negMinWP_md, PC3", 2, "cbind(log_negMinWP_md, PC3) ~ trait-1", "~us(trait):animal"),
  
  cbind("Sp_Var_covariance", "G", "log_Ks, PC1", 2, "cbind(log_Ks, PC1) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_Ks, PC2", 2, "cbind(log_Ks, PC2) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_Ks, PC3", 2, "cbind(log_Ks, PC3) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_negP50, PC1", 2, "cbind(log_negP50, PC1) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_negP50, PC2", 2, "cbind(log_negP50, PC2) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_negP50, PC3", 2, "cbind(log_negP50, PC3) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_Hv, PC1", 2, "cbind(log_Hv, PC1) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_Hv, PC2", 2, "cbind(log_Hv, PC2) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_Hv, PC3", 2, "cbind(log_Hv, PC3) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_negMinWP_md, PC1", 2, "cbind(log_negMinWP_md, PC1) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_negMinWP_md, PC2", 2, "cbind(log_negMinWP_md, PC2) ~ trait-1 + trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance", "G", "log_negMinWP_md, PC3", 2, "cbind(log_negMinWP_md, PC3) ~ trait-1 + trait:group", "~us(trait):animal"),
  
  cbind("Sp_Var_covariance_env", "GI", "log_Ks, log_negP50, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GI", "log_Ks, log_Hv, PC1, PC2, PC3", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GI", "log_Ks, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GI", "log_negP50, log_Hv, PC1, PC2, PC3", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GI", "log_negP50, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GI", "log_Hv, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group", "~us(trait):animal"),

  cbind("Sp_Var_covariance_env", "GLI", "log_Ks, log_negP50, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negP50) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group_leaf", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GLI", "log_Ks, log_Hv, PC1, PC2, PC3", 2, "cbind(log_Ks, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group_leaf", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GLI", "log_Ks, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Ks, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group_leaf", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GLI", "log_negP50, log_Hv, PC1, PC2, PC3", 2, "cbind(log_negP50, log_Hv) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group_leaf", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GLI", "log_negP50, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_negP50, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group_leaf", "~us(trait):animal"),
  cbind("Sp_Var_covariance_env", "GLI", "log_Hv, log_negMinWP_md, PC1, PC2, PC3", 2, "cbind(log_Hv, log_negMinWP_md) ~ trait-1 + (trait:PC1 + trait:PC2 + trait:PC3) * trait:group_leaf", "~us(trait):animal")
), stringsAsFactors =  T)
names(Sp_multiMdls) <- c("type", "group","variables", "N_ResVars", "fix.frml","ran.frml")
name <- paste0(as.character(Sp_multiMdls[, "type"]), "_", as.character(Sp_multiMdls[, "group"]), "_", as.character(Sp_multiMdls[, "variables"]))
name <- str_replace(name, ", ", "_")
Sp_multiMdls$name <- name
Sp_multiMdls
