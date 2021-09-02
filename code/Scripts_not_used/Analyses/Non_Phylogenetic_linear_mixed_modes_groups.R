################################# NON PHYLOGENETIC LINEAR MIXED MODELS ##################################

# Pablo Sanchez
# 02-05-2019

#http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11


### PACKAGES ###

library(ape)
library(geiger)
library(MCMCglmm)

setwd("~/Pablo/data")

# RDATA
load("RData/Non_Phylogenetic_linear_mixed_multi_models_groups.RData")


### FUNCTIONS ####

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Univariate mcmcglmm Function without filogeny
mcmcglmmNonPhyFun <- function(dta, fix.frml, ran.frml, prior, nitt = nitter, burnin = burni, thin = thinni){
  # Linear relationship between respnse and cofactor, with a random effect phylo corresponding the phylogenetic effect. The argument ginverse allows us to include a custom matrix for our random effect phylo, 
  mdl <- MCMCglmm(as.formula(fix.frml), random = as.formula(ran.frml),
                  family = "gaussian", prior = prior,
                  data = dta, nitt = nitter, burnin = burni, thin = thinni) 
  rslts <- list()
  # Model
  rslts$model <- mdl
  # Autocorrelation, effective size
  rslts$autocFix <- autocorr.diag(mdl$Sol)
  rslts$autocRan <- autocorr.diag(mdl$VCV)
  rslts$heidelFix <- heidel.diag(mdl$Sol)
  rslts$heidelRan <- heidel.diag(mdl$VCV)
  rslts$effSizeFix <- effectiveSize(mdl$VCV)
  rslts$effSizeRan <- effectiveSize(mdl$Sol)
  return(rslts)
}

# R2GLMM(m) - marginal R2GLMM (eq. 26, 29 and 30 Nakawgaway 2013)
R2m_NoPhyFun <- function(mdl){
  # alternative with crebile intervals
  n <- length(mdl$VCV[, 1])
  vmVarF<-numeric(n)
  for(i in 1:n){
    Var<-var(as.vector(mdl$Sol[i,] %*% t(mdl$X)))
    vmVarF[i]<-Var}
  R2m<-vmVarF/(vmVarF+mdl$VCV[,1]+mdl$VCV[,2])
  rslts<- list()
  rslts$R2m_mode <- posterior.mode(R2m)
  rslts$R2m_mean <- mean(R2m)
  rslts$R2m_HPD <- HPDinterval(R2m)
  return(rslts)
}

# R2GLMM(c) - conditional R2GLMM for full model (equ. 30 Nakagawa 2013)
R2c_NoPhyFun <- function(mdl){
  # alternative with crebile intervals
  n <- length(mdl$VCV[, 1])
  vmVarF<-numeric(n)
  for(i in 1:n){
    Var<-var(as.vector(mdl$Sol[i,] %*% t(mdl$X)))
    vmVarF[i]<-Var}
  R2c<-(vmVarF+mdl$VCV[,1])/(vmVarF+mdl$VCV[,1]+mdl$VCV[,2])
  rslts<- list()
  rslts$R2c_mode <- posterior.mode(R2c)
  rslts$R2c_mean <- mean(R2c)
  rslts$R2c_HPD <- HPDinterval(R2c)
  return(rslts)
}

### DATA ####

# Species datset
ht <- read.csv("HydraTRY_Pablo_2019.csv", header = T)  # Includes three principal components

# Genus tree
ht.tree <- read.tree("Phylogeny_data/Ht_tree/hydratry_genus.tre") # Genus

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% ht$genus])
plmm.data <- ht[ht$genus %in% phylo$tip.label, ]
length(plmm.data$genus)
plmm.data$family <- NULL # gives problems 
length(unique(plmm.data$genus))

# Genus random effect column
plmm.data$genus.Rand <- plmm.data$genus
plmm.data$animal <- plmm.data$genus

### PARAMETERS ####

## Prior distributions for the fixed and random effects
# G for the random effects and R for the residual variance, corresponds to a inverse-Gamma distribution with shape and scape parameter equal to 0.01, which is relatively canonical. 

# Non phylogenetic mixed models

prior1 <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 0.002)))

# Non phylogenetic models (just one random)
priorMulti2 <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2)))

# Model parameters
nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

########### HT ~ ENVIRONMENT #####################################################################################################

### log_Ks ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_env.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Ks ~ (PC1 + PC2 + PC3) * group", 
                                          ran.frml = "~ genus.Rand", prior = prior1)


## Autocorrelation and convergence
plot(log_Ks_env.nonPhyMod$model)
log_Ks_env.nonPhyMod$autocFix
log_Ks_env.nonPhyMod$autocRan
log_Ks_env.nonPhyMod$heidelFix
log_Ks_env.nonPhyMod$heidelRan
log_Ks_env.nonPhyMod$effSizeFix
log_Ks_env.nonPhyMod$effSizeRan

# Model results
summary(log_Ks_env.nonPhyMod$model)
R2c_NoPhyFun(log_Ks_env.nonPhyMod$model)
R2m_NoPhyFun(log_Ks_env.nonPhyMod$model)

### log_negP50 ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_env.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negP50 ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ genus.Rand", 
                                              prior = prior1)

## Autocorrelation and convergence
plot(log_negP50_env.nonPhyMod$model)
log_negP50_env.nonPhyMod$autocFix
log_negP50_env.nonPhyMod$autocRan
log_negP50_env.nonPhyMod$heidelFix
log_negP50_env.nonPhyMod$heidelRan
log_negP50_env.nonPhyMod$effSizeFix
log_negP50_env.nonPhyMod$effSizeRan

# Model results
summary(log_negP50_env.nonPhyMod$model)
R2c_NoPhyFun(log_negP50_env.nonPhyMod$model)
R2m_NoPhyFun(log_negP50_env.nonPhyMod$model)

### log_negMinWP_md ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_env.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negMinWP_md ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ genus.Rand", 
                                                   prior = prior1)

## Autocorrelation and convergence
plot(log_negMinWP_md_env.nonPhyMod$model)
log_negMinWP_md_env.nonPhyMod$autocFix
log_negMinWP_md_env.nonPhyMod$autocRan
log_negMinWP_md_env.nonPhyMod$heidelFix
log_negMinWP_md_env.nonPhyMod$heidelRan
log_negMinWP_md_env.nonPhyMod$effSizeFix
log_negMinWP_md_env.nonPhyMod$effSizeRan

# Model results
summary(log_negMinWP_md_env.nonPhyMod$model)
R2c_NoPhyFun(log_negMinWP_md_env.nonPhyMod$model)
R2m_NoPhyFun(log_negMinWP_md_env.nonPhyMod$model)

### log_Hv ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_env.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Hv ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ genus.Rand", 
                                          prior = prior1)

## Autocorrelation and convergence
plot(log_Hv_env.nonPhyMod$model)
log_Hv_env.nonPhyMod$autocFix
log_Hv_env.nonPhyMod$autocRan
log_Hv_env.nonPhyMod$heidelFix
log_Hv_env.nonPhyMod$heidelRan
log_Hv_env.nonPhyMod$effSizeFix
log_Hv_env.nonPhyMod$effSizeRan

# Model results
summary(log_Hv_env.nonPhyMod$model)
R2c_NoPhyFun(log_Hv_env.nonPhyMod$model)
R2m_NoPhyFun(log_Hv_env.nonPhyMod$model)

########### HT ~ HT ##############################################################################################################

### log_Ks ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_negP50.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Ks ~ log_negP50", ran.frml = "~ genus.Rand", 
                                                 prior = prior1)

## Autocorrelation and convergence
plot(log_Ks_log_negP50.nonPhyMod$model)
log_Ks_log_negP50.nonPhyMod$autocFix
log_Ks_log_negP50.nonPhyMod$autocRan
log_Ks_log_negP50.nonPhyMod$heidelFix
log_Ks_log_negP50.nonPhyMod$heidelRan
log_Ks_log_negP50.nonPhyMod$effSizeFix
log_Ks_log_negP50.nonPhyMod$effSizeRan

# Model results
summary(log_Ks_log_negP50.nonPhyMod$model)

r2m <- R2m_NoPhyFun(log_Ks_log_negP50.nonPhyMod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2c_NoPhyFun(log_Ks_log_negP50.nonPhyMod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negP50 ~ log_Ks ###

# Model
log_negP50_log_Ks.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negP50 ~ log_Ks", ran.frml = "~ genus.Rand", 
                                                 prior = prior1)

## Autocorrelation and convergence
plot(log_negP50_log_Ks.nonPhyMod$model)
log_negP50_log_Ks.nonPhyMod$autocFix
log_negP50_log_Ks.nonPhyMod$autocRan
log_negP50_log_Ks.nonPhyMod$heidelFix
log_negP50_log_Ks.nonPhyMod$heidelRan
log_negP50_log_Ks.nonPhyMod$effSizeFix
log_negP50_log_Ks.nonPhyMod$effSizeRan

# Model results
summary(log_negP50_log_Ks.nonPhyMod$model)

r2m <- R2m_NoPhyFun(log_negP50_log_Ks.nonPhyMod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2c_NoPhyFun(log_negP50_log_Ks.nonPhyMod$model)
r2c$R2c_mean
r2c$R2c_HPD


### log_negMinWP_md ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negMinWP_md_log_negP50.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negMinWP_md ~ log_negP50", ran.frml = "~ genus.Rand", 
                                                          prior = prior1)

## Autocorrelation and convergence
plot(log_negMinWP_md_log_negP50.nonPhyMod$model)
log_negMinWP_md_log_negP50.nonPhyMod$autocFix
log_negMinWP_md_log_negP50.nonPhyMod$autocRan
log_negMinWP_md_log_negP50.nonPhyMod$heidelFix
log_negMinWP_md_log_negP50.nonPhyMod$heidelRan
log_negMinWP_md_log_negP50.nonPhyMod$effSizeFix
log_negMinWP_md_log_negP50.nonPhyMod$effSizeRan

# Model results
summary(log_negMinWP_md_log_negP50.nonPhyMod$model)

r2m <- R2m_NoPhyFun(log_negMinWP_md_log_negP50.nonPhyMod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2c_NoPhyFun(log_negMinWP_md_log_negP50.nonPhyMod$model)
r2c$R2c_mean
r2c$R2c_HPD


### log_negP50 ~ log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negP50_log_negMinWP_md.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negP50 ~ log_negMinWP_md * group", ran.frml = "~ genus.Rand", 
                                                          prior = prior1)

## Autocorrelation and convergence
plot(log_negP50_log_negMinWP_md.nonPhyMod$model)
log_negP50_log_negMinWP_md.nonPhyMod$autocFix
log_negP50_log_negMinWP_md.nonPhyMod$autocRan
log_negP50_log_negMinWP_md.nonPhyMod$heidelFix
log_negP50_log_negMinWP_md.nonPhyMod$heidelRan
log_negP50_log_negMinWP_md.nonPhyMod$effSizeFix
log_negP50_log_negMinWP_md.nonPhyMod$effSizeRan

# Model results
summary(log_negP50_log_negMinWP_md.nonPhyMod$model)

r2m <- R2m_NoPhyFun(log_negMinWP_md_log_negP50.nonPhyMod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2c_NoPhyFun(log_negMinWP_md_log_negP50.nonPhyMod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Hv ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_log_Ks.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Hv ~ log_Ks", ran.frml = "~ genus.Rand", 
                                             prior = prior1)

## Autocorrelation and convergence
plot(log_Hv_log_Ks.nonPhyMod$model)
log_Hv_log_Ks.nonPhyMod$autocFix
log_Hv_log_Ks.nonPhyMod$autocRan
log_Hv_log_Ks.nonPhyMod$heidelFix
log_Hv_log_Ks.nonPhyMod$heidelRan
log_Hv_log_Ks.nonPhyMod$effSizeFix
log_Hv_log_Ks.nonPhyMod$effSizeRan

# Model results
summary(log_Hv_log_Ks.nonPhyMod$model)

r2m <- R2m_NoPhyFun(log_Hv_log_Ks.nonPhyMod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2c_NoPhyFun(log_Hv_log_Ks.nonPhyMod$model)
r2c$R2c_mean
r2c$R2c_HPD


### log_Ks ~ log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_Hv.nonPhyMod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Ks ~ log_Hv", ran.frml = "~ genus.Rand", 
                                             prior = prior1)

## Autocorrelation and convergence
plot(log_Ks_log_Hv.nonPhyMod$model)
log_Ks_log_Hv.nonPhyMod$autocFix
log_Ks_log_Hv.nonPhyMod$autocRan
log_Ks_log_Hv.nonPhyMod$heidelFix
log_Ks_log_Hv.nonPhyMod$heidelRan
log_Ks_log_Hv.nonPhyMod$effSizeFix
log_Ks_log_Hv.nonPhyMod$effSizeRan

# Model results
summary(log_Ks_log_Hv.nonPhyMod$model)

r2m <- R2m_NoPhyFun(log_Ks_log_Hv.nonPhyMod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2c_NoPhyFun(log_Ks_log_Hv.nonPhyMod$model)
r2c$R2c_mean
r2c$R2c_HPD


########### HT ~ HT * group ##############################################################################################################


### log_Ks ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_negP50_gr_NoPhy.mod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Ks ~ log_negP50 * group", ran.frml = "~ genus.Rand", 
                                                 prior = prior1)

## Autocorrelation and convergence
# plot(log_Ks_log_negP50_NoPhy.mod$model)
log_Ks_log_negP50_NoPhy.mod$autocFix
log_Ks_log_negP50_NoPhy.mod$autocRan
log_Ks_log_negP50_NoPhy.mod$heidelFix
log_Ks_log_negP50_NoPhy.mod$heidelRan
log_Ks_log_negP50_NoPhy.mod$effSizeFix
log_Ks_log_negP50_NoPhy.mod$effSizeRan

# Model results
summary(log_Ks_log_negP50_NoPhy.mod$model)

### log_negP50 ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_log_Ks_gr_NoPhy.mod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negP50 ~ log_Ks * group", ran.frml = "~ genus.Rand", 
                                                 prior = prior1)

## Autocorrelation and convergence
# plot(log_negP50_log_Ks_NoPhy.mod$model)
log_negP50_log_Ks_NoPhy.mod$autocFix
log_negP50_log_Ks_NoPhy.mod$autocRan
log_negP50_log_Ks_NoPhy.mod$heidelFix
log_negP50_log_Ks_NoPhy.mod$heidelRan
log_negP50_log_Ks_NoPhy.mod$effSizeFix
log_negP50_log_Ks_NoPhy.mod$effSizeRan

# Model results
summary(log_negP50_log_Ks_NoPhy.mod$model)


### log_negP50 ~ log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negP50_log_negMinWP_md_gr_NoPhy.mod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negP50 ~ log_negMinWP_md * group", ran.frml = "~ genus.Rand", 
                                                          prior = prior1)

## Autocorrelation and convergence
# plot(log_negP50_log_negMinWP_md_NoPhy.mod$model)
log_negP50_log_negMinWP_md_NoPhy.mod$autocFix
log_negP50_log_negMinWP_md_NoPhy.mod$autocRan
log_negP50_log_negMinWP_md_NoPhy.mod$heidelFix
log_negP50_log_negMinWP_md_NoPhy.mod$heidelRan
log_negP50_log_negMinWP_md_NoPhy.mod$effSizeFix
log_negP50_log_negMinWP_md_NoPhy.mod$effSizeRan

# Model results
summary(log_negP50_log_negMinWP_md_NoPhy.mod$model)


### log_negMinWP_md ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negMinWP_md_log_negP50_gr_NoPhy.mod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_negMinWP_md ~ log_negP50", ran.frml = "~ genus.Rand", 
                                                          prior = prior1)


## Autocorrelation and convergence
# plot(log_negMinWP_md_log_negP50_NoPhy.mod.mod$model)
log_negMinWP_md_log_negP50_NoPhy.mod.mod$autocFix
log_negMinWP_md_log_negP50_NoPhy.mod.mod$autocRan
log_negMinWP_md_log_negP50_NoPhy.mod.mod$heidelFix
log_negMinWP_md_log_negP50_NoPhy.mod.mod$heidelRan
log_negMinWP_md_log_negP50_NoPhy.mod.mod$effSizeFix
log_negMinWP_md_log_negP50_NoPhy.mod.mod$effSizeRan

# Model results
summary(log_negMinWP_md_log_negP50_NoPhy.mod$model)

### log_Hv ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_log_Ks_gr_NoPhy.mod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Hv ~ log_Ks", ran.frml = "~ genus.Rand", 
                                             prior = prior1)

## Autocorrelation and convergence
# plot(log_Hv_log_Ks_NoPhy.mod$model)
log_Hv_log_Ks_NoPhy.mod$autocFix
log_Hv_log_Ks_NoPhy.mod$autocRan
log_Hv_log_Ks_NoPhy.mod$heidelFix
log_Hv_log_Ks_NoPhy.mod$heidelRan
log_Hv_log_Ks_NoPhy.mod$effSizeFix
log_Hv_log_Ks_NoPhy.mod$effSizeRan


# Model results
summary(log_Hv_log_Ks_NoPhy.mod$model)


### log_Ks ~ log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_Hv_gr_NoPhy.mod <- mcmcglmmNonPhyFun(dta = data, fix.frml = "log_Ks ~ log_Hv", ran.frml = "~ genus.Rand", 
                                             prior = prior1)

## Autocorrelation and convergence
# plot(log_Ks_log_Hv_NoPhy.mod$model)
log_Ks_log_Hv_NoPhy.mod$autocFix
log_Ks_log_Hv_NoPhy.mod$autocRan
log_Ks_log_Hv_NoPhy.mod$heidelFix
log_Ks_log_Hv_NoPhy.mod$heidelRan
log_Ks_log_Hv_NoPhy.mod$effSizeFix
log_Ks_log_Hv_NoPhy.mod$effSizeRan


# Model results
summary(log_Ks_log_Hv_NoPhy.mod$model)


############################# HT + HT ~ -1 ##############################################################################

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
meanKsP50.nonPhyMultiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:group, 
                                     random = ~us(trait):genus.Rand, 
                                     rcov = ~us(trait):units, family = c("gaussian","gaussian"),
                                     prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsP50.nonPhyMultiMod)
plot(meanKsP50.nonPhyMultiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsP50.nonPhyMultiMod$Sol)
autocorr.diag(meanKsP50.nonPhyMultiMod$VCV)
# Convergence
heidel.diag(meanKsP50.nonPhyMultiMod$Sol)
heidel.diag(meanKsP50.nonPhyMultiMod$VCV)

## Evolutive correlation
corr.evo <- meanKsP50.nonPhyMultiMod$VCV[, "traitlog_negP50:traitlog_Ks.units"]/
  sqrt(meanKsP50.nonPhyMultiMod$VCV[, "traitlog_Ks:traitlog_Ks.units"]*meanKsP50.nonPhyMultiMod$VCV[,"traitlog_negP50:traitlog_negP50.units"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
meanKsHv.nonPhyMultiMod <- MCMCglmm(cbind(log_Ks, log_Hv) ~ trait-1  + trait:group, 
                                    random = ~ us(trait):genus.Rand, 
                                    rcov = ~us(trait):units, family = c("gaussian","gaussian"),
                                    prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsHv.nonPhyMultiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsHv.nonPhyMultiMod$Sol)
autocorr.diag(meanKsHv.nonPhyMultiMod$VCV)
# Convergence
heidel.diag(meanKsHv.nonPhyMultiMod$Sol)
heidel.diag(meanKsHv.nonPhyMultiMod$VCV)

# Evolutive correlation
corr.evo <- meanKsHv.nonPhyMultiMod$VCV[, "traitlog_Hv:traitlog_Ks.units"]/
  sqrt(meanKsHv.nonPhyMultiMod$VCV[, "traitlog_Ks:traitlog_Ks.units"]*meanKsHv.nonPhyMultiMod$VCV[,"traitlog_Hv:traitlog_Hv.units"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
meanMinWPP50.nonPhyMultiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1  + trait:group, 
                                        random = ~ us(trait):genus.Rand, 
                                        rcov = ~us(trait):units, family = c("gaussian","gaussian"),
                                        prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanMinWPP50.nonPhyMultiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanMinWPP50.nonPhyMultiMod$Sol)
autocorr.diag(meanMinWPP50.nonPhyMultiMod$VCV)
# Convergence
heidel.diag(meanMinWPP50.nonPhyMultiMod$Sol)
heidel.diag(meanMinWPP50.nonPhyMultiMod$VCV)

# Evolutive correlation
corr.evo <- meanMinWPP50.nonPhyMultiMod$VCV[, "traitlog_negMinWP_md:traitlog_negP50.units"]/
  sqrt(meanMinWPP50.nonPhyMultiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.units"]*meanMinWPP50.nonPhyMultiMod$VCV[,"traitlog_negP50:traitlog_negP50.units"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)


############################# HT + HT ~ ENVIRONMENT ##############################################################################

# Complete data
plmm.data <- completeFun(plmm.data, c("PC1", "Species"))

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50.nonPhyMultiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                             random = ~us(trait):genus.Rand, 
                                             rcov = ~us(trait):units, family = c("gaussian","gaussian"),
                                             prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50.nonPhyMultiMod)
plot(log_Ks_log_negP50.nonPhyMultiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50.nonPhyMultiMod$Sol)
autocorr.diag(log_Ks_log_negP50.nonPhyMultiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50.nonPhyMultiMod$Sol)
heidel.diag(log_Ks_log_negP50.nonPhyMultiMod$VCV)


## Evolutive correlation
corr.evo <- log_Ks_log_negP50.nonPhyMultiMod$VCV[, "traitlog_negP50:traitlog_Ks.units"]/
  sqrt(log_Ks_log_negP50.nonPhyMultiMod$VCV[, "traitlog_Ks:traitlog_Ks.units"]*log_Ks_log_negP50.nonPhyMultiMod$VCV[,"traitlog_negP50:traitlog_negP50.units"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_Hv.nonPhyMultiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                         random = ~ us(trait):genus.Rand, 
                                         rcov = ~us(trait):units, family = c("gaussian","gaussian"),
                                         prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv.nonPhyMultiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv.nonPhyMultiMod$Sol)
autocorr.diag(log_Ks_log_Hv.nonPhyMultiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv.nonPhyMultiMod$Sol)
heidel.diag(log_Ks_log_Hv.nonPhyMultiMod$VCV)


# Evolutive correlation
corr.evo <- log_Ks_log_Hv.nonPhyMultiMod$VCV[, "traitlog_Hv:traitlog_Ks.units"]/
  sqrt(log_Ks_log_Hv.nonPhyMultiMod$VCV[, "traitlog_Ks:traitlog_Ks.units"]*log_Ks_log_Hv.nonPhyMultiMod$VCV[,"traitlog_Hv:traitlog_Hv.units"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)

### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_negP50.nonPhyMultiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                                      random = ~ us(trait):genus.Rand, 
                                                      rcov = ~us(trait):units, family = c("gaussian","gaussian"),
                                                      prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_negP50.nonPhyMultiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_negP50.nonPhyMultiMod$Sol)
autocorr.diag(log_negMinWP_md_log_negP50.nonPhyMultiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_negP50.nonPhyMultiMod$Sol)
heidel.diag(log_negMinWP_md_log_negP50.nonPhyMultiMod$VCV)


# Evolutive correlation
corr.evo <- log_negMinWP_md_log_negP50.nonPhyMultiMod$VCV[, "traitlog_negMinWP_md:traitlog_negP50.units"]/
  sqrt(log_negMinWP_md_log_negP50.nonPhyMultiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.units"]*log_negMinWP_md_log_negP50.nonPhyMultiMod$VCV[,"traitlog_negP50:traitlog_negP50.units"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)


#### SAVE IMAGE ####
save.image(file = "RData/Non_Phylogenetic_linear_mixed_multi_models_groups.RData")
