################################# PHYLOGENETIC LINEAR MIXED MODELS ##################################

# Pablo Sanchez
# 02-05-2019

#http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11

### PACKAGES ###
# install.packages("rstanarm")

library(ape)
library(geiger)
library(MCMCglmm)

setwd("~/Pablo/data")

# RDATA
load("RData/not_used/Phylogenetic_linear_mixed_multi_models.RData")


### FUNCTIONS ####

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Univariate mcmcglmm Function
mcmcglmmFun <- function(dta, fix.frml, ran.frml, prior, nitte = nitter, burn = burni, thinn = thinni){
  # Linear relationship between respnse and cofactor, with a random effect phylo corresponding the phylogenetic effect. The argument ginverse allows us to include a custom matrix for our random effect phylo, 
  mdl <- MCMCglmm(as.formula(fix.frml), random = as.formula(ran.frml),
                  family = "gaussian", pedigree=phylo, prior = prior,
                  data = dta, nitt = nitte, burnin = burn, thin = thinn) 
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
  ## Phylogenetic signal
  lambda <- mdl$VCV[,'animal']/
    (mdl$VCV[,'animal'] + mdl$VCV[,'genus.Rand'] + mdl$VCV[,'units'])
  rslts$lambda <- lambda
  rslts$lambdaMean <- mean(lambda)
  rslts$lambdaPost.mode <- posterior.mode(lambda)
  rslts$lambdaHPDinterval <- HPDinterval(lambda) # 95% credible interval 
  return(rslts)
}


# R2GLMM(m) - marginal R2GLMM (eq. 26, 29 and 30 Nakawgaway 2013)
R2mFun <- function(mdl){
  # alternative with crebile intervals
  n <- length(mdl$VCV[, 1])
  vmVarF<-numeric(n)
  for(i in 1:n){
    Var<-var(as.vector(mdl$Sol[i,] %*% t(mdl$X)))
    vmVarF[i]<-Var}
  R2m<-vmVarF/(vmVarF+mdl$VCV[,1]+mdl$VCV[,2]+mdl$VCV[,3])
  rslts<- list()
  rslts$R2m_mode <- posterior.mode(R2m)
  rslts$R2m_mean <- mean(R2m)
  rslts$R2m_HPD <- HPDinterval(R2m)
  return(rslts)
}

# R2GLMM(c) - conditional R2GLMM for full model (equ. 30 Nakagawa 2013)
R2cFun <- function(mdl){
  # alternative with crebile intervals
  n <- length(mdl$VCV[, 1])
  vmVarF<-numeric(n)
  for(i in 1:n){
    Var<-var(as.vector(mdl$Sol[i,] %*% t(mdl$X)))
    vmVarF[i]<-Var}
  R2c<-(vmVarF+mdl$VCV[,1]+mdl$VCV[,2])/(vmVarF+mdl$VCV[,1]+mdl$VCV[,2]+mdl$VCV[,3])
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

# Phylogenetic mixed models
prior2 <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 0.002), 
           G2 = list(V = 1, nu = 0.002))) # two random effects, phylogeny and genus

# Phylogenetic mixed multi models
priorMulti <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2), 
         G2 = list(V=diag(2)/2,nu=2)))

# Model parameters
nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

########### PHYLOGENETIC SIGNAL #####################################################################################################

### PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
PC1.mod <- mcmcglmmFun(dta = data, fix.frml = "PC1 ~ 1", ran.frml = "~ animal + genus.Rand", 
                       prior = prior2)

## Autocorrelation and convergence
plot(PC1.mod$model)
PC1.mod$autocFix
PC1.mod$autocRan
PC1.mod$heidelFix
PC1.mod$heidelRan
PC1.mod$effSizeFix
PC1.mod$effSizeRan

# Lambda
plot(PC1.mod$lambda)
PC1.mod$lambdaMean
PC1.mod$lambdaHPDinterval

lambda <- PC1.mod$model$VCV[,'animal']/
  (PC1.mod$model$VCV[,'animal'] + PC1.mod$model$VCV[,'genus.Rand'] + PC1.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Genus lability
lab.PC1 <- PC1.mod$model$VCV[,'genus.Rand']/
  (PC1.mod$model$VCV[,'animal'] + PC1.mod$model$VCV[,'genus.Rand'] + PC1.mod$model$VCV[,'units'])
mean(lab.PC1)
HPDinterval(lab.PC1, 0.95)

# Intra-genus variation
intra.PC1 <- PC1.mod$model$VCV[,'units']/
  (PC1.mod$model$VCV[,'animal'] + PC1.mod$model$VCV[,'genus.Rand'] + PC1.mod$model$VCV[,'units'])
mean(intra.PC1)
HPDinterval(intra.PC1, 0.95)

# Model results
summary(PC1.mod$model)
R2cFun(PC1.mod$model)

### PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
PC2.mod <- mcmcglmmFun(dta = data, fix.frml = "PC2 ~ 1", ran.frml = "~ animal + genus.Rand", 
                       prior = prior2)

## Autocorrelation and convergence
plot(PC2.mod$model)
PC2.mod$autocFix
PC2.mod$autocRan
PC2.mod$heidelFix
PC2.mod$heidelRan
PC2.mod$effSizeFix
PC2.mod$effSizeRan

# Lambda
plot(PC2.mod$lambda)
PC2.mod$lambdaMean
PC2.mod$lambdaHPDinterval

lambda <- PC2.mod$model$VCV[,'animal']/
  (PC2.mod$model$VCV[,'animal'] + PC2.mod$model$VCV[,'genus.Rand'] + PC2.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability

lab.PC2 <- PC2.mod$model$VCV[,'genus.Rand']/
  (PC2.mod$model$VCV[,'animal'] + PC2.mod$model$VCV[,'genus.Rand'] + PC2.mod$model$VCV[,'units'])
mean(lab.PC2)
HPDinterval(lab.PC2, 0.95)

# Intra-genus

intra.PC2 <- PC2.mod$model$VCV[,'units']/
  (PC2.mod$model$VCV[,'animal'] + PC2.mod$model$VCV[,'genus.Rand'] + PC2.mod$model$VCV[,'units'])
mean(intra.PC2)
HPDinterval(intra.PC2, 0.95)

# Model results
summary(PC2.mod$model)
R2cFun(PC2.mod$model)

### PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
PC3.mod <- mcmcglmmFun(dta = data, fix.frml = "PC3 ~ 1", ran.frml = "~ animal + genus.Rand", 
                       prior = prior2)

## Autocorrelation and convergence
plot(PC3.mod$model)
PC3.mod$autocFix
PC3.mod$autocRan
PC3.mod$heidelFix
PC3.mod$heidelRan
PC3.mod$effSizeFix
PC3.mod$effSizeRan

# Lambda
plot(PC3.mod$lambda)
PC3.mod$lambdaMean
PC3.mod$lambdaHPDinterval

lambda <- PC3.mod$model$VCV[,'animal']/
  (PC3.mod$model$VCV[,'animal'] + PC3.mod$model$VCV[,'genus.Rand'] + PC3.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# lability

lab.PC3 <- PC3.mod$model$VCV[,'genus.Rand']/
  (PC3.mod$model$VCV[,'animal'] + PC3.mod$model$VCV[,'genus.Rand'] + PC3.mod$model$VCV[,'units'])
mean(lab.PC3)
HPDinterval(lab.PC3, 0.95)

# Intra-genus

intra.PC3 <- PC3.mod$model$VCV[,'units']/
  (PC3.mod$model$VCV[,'animal'] + PC3.mod$model$VCV[,'genus.Rand'] + PC3.mod$model$VCV[,'units'])
mean(intra.PC3)
HPDinterval(intra.PC3, 0.95)

# Model results
summary(PC3.mod$model)
R2cFun(PC3.mod$model)

### log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(log_Ks.mod$model)
log_Ks.mod$autocFix
log_Ks.mod$autocRan
log_Ks.mod$heidelFix
log_Ks.mod$heidelRan
log_Ks.mod$effSizeFix
log_Ks.mod$effSizeRan

# Lambda
plot(log_Ks.mod$lambda)
log_Ks.mod$lambdaMean
log_Ks.mod$lambdaHPDinterval

lambda <- log_Ks.mod$model$VCV[,'animal']/
  (log_Ks.mod$model$VCV[,'animal'] + log_Ks.mod$model$VCV[,'genus.Rand'] + log_Ks.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# lability

lab.Ks <- log_Ks.mod$model$VCV[,'genus.Rand']/
  (log_Ks.mod$model$VCV[,'animal'] + log_Ks.mod$model$VCV[,'genus.Rand'] + log_Ks.mod$model$VCV[,'units'])
mean(lab.Ks)
HPDinterval(lab.Ks, 0.95)

# Intra-genus

intra.Ks <- log_Ks.mod$model$VCV[,'units']/
  (log_Ks.mod$model$VCV[,'animal'] + log_Ks.mod$model$VCV[,'genus.Rand'] + log_Ks.mod$model$VCV[,'units'])
mean(intra.Ks)
HPDinterval(intra.Ks, 0.95)

# Model results
summary(log_Ks.mod$model)
R2cFun(log_Ks.mod$model)

### log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ 1", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_negP50.mod$model)
log_negP50.mod$autocFix
log_negP50.mod$autocRan
log_negP50.mod$heidelFix
log_negP50.mod$heidelRan
log_negP50.mod$effSizeFix
log_negP50.mod$effSizeRan

# Lambda
plot(log_negP50.mod$lambda)
log_negP50.mod$lambdaMean
log_negP50.mod$lambdaHPDinterval

lambda <- log_negP50.mod$model$VCV[,'animal']/
  (log_negP50.mod$model$VCV[,'animal'] + log_negP50.mod$model$VCV[,'genus.Rand'] + log_negP50.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability

lab.P50 <- log_negP50.mod$model$VCV[,'genus.Rand']/
  (log_negP50.mod$model$VCV[,'animal'] + log_negP50.mod$model$VCV[,'genus.Rand'] + log_negP50.mod$model$VCV[,'units'])
mean(lab.P50)
HPDinterval(lab.P50, 0.95)

# Intra-genus

intra.P50 <- log_negP50.mod$model$VCV[,'units']/
  (log_negP50.mod$model$VCV[,'animal'] + log_negP50.mod$model$VCV[,'genus.Rand'] + log_negP50.mod$model$VCV[,'units'])
mean(intra.P50)
HPDinterval(intra.P50, 0.95)

# Model results
summary(log_negP50.mod$model)
R2cFun(log_negP50.mod$model)

### log_minWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negMinWP_md", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ 1", ran.frml = "~ animal + genus.Rand", 
                                   prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md.mod$model)
log_negMinWP_md.mod$autocFix
log_negMinWP_md.mod$autocRan
log_negMinWP_md.mod$heidelFix
log_negMinWP_md.mod$heidelRan
log_negMinWP_md.mod$effSizeFix
log_negMinWP_md.mod$effSizeRan

# Lambda
plot(log_negMinWP_md.mod$lambda)
log_negMinWP_md.mod$lambdaMean
log_negMinWP_md.mod$lambdaHPDinterval

lambda <- log_negMinWP_md.mod$model$VCV[,'animal']/
  (log_negMinWP_md.mod$model$VCV[,'animal'] + log_negMinWP_md.mod$model$VCV[,'genus.Rand'] + log_negMinWP_md.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability

lab.MinWP <- log_negMinWP_md.mod$model$VCV[,'genus.Rand']/
  (log_negMinWP_md.mod$model$VCV[,'animal'] + log_negMinWP_md.mod$model$VCV[,'genus.Rand'] + log_negMinWP_md.mod$model$VCV[,'units'])
mean(lab.MinWP)
HPDinterval(lab.MinWP, 0.95)

# Intra-genus 

intra.MinWP <- log_negMinWP_md.mod$model$VCV[,'units']/
  (log_negMinWP_md.mod$model$VCV[,'animal'] + log_negMinWP_md.mod$model$VCV[,'genus.Rand'] + log_negMinWP_md.mod$model$VCV[,'units'])
mean(intra.MinWP)
HPDinterval(intra.MinWP, 0.95)

# Model results
summary(log_negMinWP_md.mod$model)
R2cFun(log_negMinWP_md.mod$model)

### log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(log_Hv.mod$model)
log_Hv.mod$autocFix
log_Hv.mod$autocRan
log_Hv.mod$heidelFix
log_Hv.mod$heidelRan
log_Hv.mod$effSizeFix
log_Hv.mod$effSizeRan

# Lambda
plot(log_Hv.mod$lambda)
log_Hv.mod$lambdaMean
log_Hv.mod$lambdaHPDinterval

lambda <- log_Hv.mod$model$VCV[,'animal']/
  (log_Hv.mod$model$VCV[,'animal'] + log_Hv.mod$model$VCV[,'genus.Rand'] + log_Hv.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability

lab.Hv <- log_Hv.mod$model$VCV[,'genus.Rand']/
  (log_Hv.mod$model$VCV[,'animal'] + log_Hv.mod$model$VCV[,'genus.Rand'] + log_Hv.mod$model$VCV[,'units'])
mean(lab.Hv)
HPDinterval(lab.Hv, 0.95)

# Intra-genus

intra.Hv <- log_Hv.mod$model$VCV[,'units']/
  (log_Hv.mod$model$VCV[,'animal'] + log_Hv.mod$model$VCV[,'genus.Rand'] + log_Hv.mod$model$VCV[,'units'])
mean(intra.Hv)
HPDinterval(intra.Hv, 0.95)

# Model results
summary(log_Hv.mod$model)
R2cFun(log_Hv.mod$model)

### log_negPItlp ###

# Complete data
data <- completeFun(plmm.data, c("log_negPItlp", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negPItlp.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negPItlp ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(log_negPItlp.mod$model)
log_negPItlp.mod$autocFix
log_negPItlp.mod$autocRan
log_negPItlp.mod$heidelFix
log_negPItlp.mod$heidelRan
log_negPItlp.mod$effSizeFix
log_negPItlp.mod$effSizeRan

# Lambda
plot(log_negPItlp.mod$lambda)
log_negPItlp.mod$lambdaMean
log_negPItlp.mod$lambdaHPDinterval

lambda <- log_negPItlp.mod$model$VCV[,'animal']/
  (log_negPItlp.mod$model$VCV[,'animal'] + log_negPItlp.mod$model$VCV[,'genus.Rand'] + log_negPItlp.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability 

lab.Ptlp <- log_negPItlp.mod$model$VCV[,'genus.Rand']/
  (log_negPItlp.mod$model$VCV[,'animal'] + log_negPItlp.mod$model$VCV[,'genus.Rand'] + log_negPItlp.mod$model$VCV[,'units'])
mean(lab.Ptlp)
HPDinterval(lab.Ptlp, 0.95)

# Intra-genus

intra.Ptlp <- log_negPItlp.mod$model$VCV[,'units']/
  (log_negPItlp.mod$model$VCV[,'animal'] + log_negPItlp.mod$model$VCV[,'genus.Rand'] + log_negPItlp.mod$model$VCV[,'units'])
mean(intra.Ptlp)
HPDinterval(intra.Ptlp, 0.95)

# Model results
summary(log_negPItlp.mod$model)


### log_Kl ###

# Complete data
data <- completeFun(plmm.data, c("log_Kl", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Kl.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Kl ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(log_Kl.mod$model)
log_Kl.mod$autocFix
log_Kl.mod$autocRan
log_Kl.mod$heidelFix
log_Kl.mod$heidelRan
log_Kl.mod$effSizeFix
log_Kl.mod$effSizeRan

# Lambda
plot(log_Kl.mod$lambda)
log_Kl.mod$lambdaMean
log_Kl.mod$lambdaHPDinterval

lambda <- log_Kl.mod$model$VCV[,'animal']/
  (log_Kl.mod$model$VCV[,'animal'] + log_Kl.mod$model$VCV[,'genus.Rand'] + log_Kl.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability

lab.Kl <- log_Kl.mod$model$VCV[,'genus.Rand']/
  (log_Kl.mod$model$VCV[,'animal'] + log_Kl.mod$model$VCV[,'genus.Rand'] + log_Kl.mod$model$VCV[,'units'])
mean(lab.Kl)
HPDinterval(lab.Kl, 0.95)

# Intra-genus

intra.Kl <- log_Kl.mod$model$VCV[,'units']/
  (log_Kl.mod$model$VCV[,'animal'] + log_Kl.mod$model$VCV[,'genus.Rand'] + log_Kl.mod$model$VCV[,'units'])
mean(intra.Kl)
HPDinterval(intra.Kl, 0.95)

# Model results
summary(log_Kl.mod$model)

### log_Hmax ###

# Complete data
data <- completeFun(plmm.data, c("log_Hmax", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hmax.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hmax ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(log_Hmax.mod$model)
log_Hmax.mod$autocFix
log_Hmax.mod$autocRan
log_Hmax.mod$heidelFix
log_Hmax.mod$heidelRan
log_Hmax.mod$effSizeFix
log_Hmax.mod$effSizeRan

# Lambda
plot(log_Hmax.mod$lambda)
log_Hmax.mod$lambdaMean
log_Hmax.mod$lambdaHPDinterval

lambda <- log_Hmax.mod$model$VCV[,'animal']/
  (log_Hmax.mod$model$VCV[,'animal'] + log_Hmax.mod$model$VCV[,'genus.Rand'] + log_Hmax.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Lability

lab.Hmax <- log_Hmax.mod$model$VCV[,'genus.Rand']/
  (log_Hmax.mod$model$VCV[,'animal'] + log_Hmax.mod$model$VCV[,'genus.Rand'] + log_Hmax.mod$model$VCV[,'units'])
mean(lab.Hmax)
HPDinterval(lab.Hmax, 0.95)


# Model results
summary(log_Hmax.mod$model)

### Rd ###

# Complete data
data <- completeFun(plmm.data, c("Rd", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
Rd.mod <- mcmcglmmFun(dta = data, fix.frml = "Rd ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(Rd.mod$model)
Rd.mod$autocFix
Rd.mod$autocRan
Rd.mod$heidelFix
Rd.mod$heidelRan
Rd.mod$effSizeFix
Rd.mod$effSizeRan

# Lambda
plot(Rd.mod$lambda)
Rd.mod$lambdaMean
Rd.mod$lambdaHPDinterval

lambda <- Rd.mod$model$VCV[,'animal']/
  (Rd.mod$model$VCV[,'animal'] + Rd.mod$model$VCV[,'genus.Rand'] + Rd.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Model results
summary(Rd.mod$model)

### log_gmax ###

# Complete data
data <- completeFun(plmm.data, c("log_gmax", "genus"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_gmax.mod <- mcmcglmmFun(dta = data, fix.frml = "log_gmax ~ 1", ran.frml = "~ animal + genus.Rand", 
                          prior = prior2)

## Autocorrelation and convergence
plot(log_gmax.mod$model)
log_gmax.mod$autocFix
log_gmax.mod$autocRan
log_gmax.mod$heidelFix
log_gmax.mod$heidelRan
log_gmax.mod$effSizeFix
log_gmax.mod$effSizeRan

# Lambda
plot(log_gmax.mod$lambda)
log_gmax.mod$lambdaMean
log_gmax.mod$lambdaHPDinterval

lambda <- log_gmax.mod$model$VCV[,'animal']/
  (log_gmax.mod$model$VCV[,'animal'] + log_gmax.mod$model$VCV[,'genus.Rand'] + log_gmax.mod$model$VCV[,'units'])
mean(lambda)
HPDinterval(lambda, 0.95)

# Model results
summary(log_gmax.mod$model)

### PS RESULTS SUMMARY ###

PC1.ps <- data.frame("Variable" = "PC1", "lambda"  = PC1.mod$lambdaMean, 
                     "lower HPD" = PC1.mod$lambdaHPDinterval[1], 
                     "upper HPD" = PC1.mod$lambdaHPDinterval[2])
PC2.ps <- data.frame("Variable" = "PC2", "lambda"  = PC2.mod$lambdaMean, 
                     "lower HPD" = PC2.mod$lambdaHPDinterval[1], 
                     "upper HPD" = PC2.mod$lambdaHPDinterval[2])
PC3.ps <- data.frame("Variable" = "PC3", "lambda"  = PC3.mod$lambdaMean, 
                     "lower HPD" = PC3.mod$lambdaHPDinterval[1], 
                     "upper HPD" = PC3.mod$lambdaHPDinterval[2])
log_Ks.ps <- data.frame("Variable" = "log_Ks", "lambda"  = log_Ks.mod$lambdaMean, 
                        "lower HPD" = log_Ks.mod$lambdaHPDinterval[1], 
                        "upper HPD" = log_Ks.mod$lambdaHPDinterval[2])
log_Hv.ps <- data.frame("Variable" = "log_Hv", "lambda"  = log_Hv.mod$lambdaMean, 
                        "lower HPD" = log_Hv.mod$lambdaHPDinterval[1], 
                        "upper HPD" = log_Hv.mod$lambdaHPDinterval[2])
log_negP50.ps <- data.frame("Variable" = "log_negP50", "lambda"  = log_negP50.mod$lambdaMean, 
                            "lower HPD" = log_negP50.mod$lambdaHPDinterval[1], 
                            "upper HPD" = log_negP50.mod$lambdaHPDinterval[2])
log_negMinWP_md.ps <- data.frame("Variable" = "log_negMinWP_md", "lambda"  = log_negMinWP_md.mod$lambdaMean, 
                                 "lower HPD" = log_negMinWP_md.mod$lambdaHPDinterval[1], 
                                 "upper HPD" = log_negMinWP_md.mod$lambdaHPDinterval[2])
log_negPItlp.ps <- data.frame("Variable" = "log_negPItlp", "lambda"  = log_negPItlp.mod$lambdaMean, 
                        "lower HPD" = log_negPItlp.mod$lambdaHPDinterval[1], 
                        "upper HPD" = log_negPItlp.mod$lambdaHPDinterval[2])
log_Kl.ps <- data.frame("Variable" = "log_Kl", "lambda"  = log_Kl.mod$lambdaMean, 
                        "lower HPD" = log_Kl.mod$lambdaHPDinterval[1], 
                        "upper HPD" = log_Kl.mod$lambdaHPDinterval[2])
log_Hmax.ps <- data.frame("Variable" = "log_Hmax", "lambda"  = log_Hmax.mod$lambdaMean, 
                        "lower HPD" = log_Hmax.mod$lambdaHPDinterval[1], 
                        "upper HPD" = log_Hmax.mod$lambdaHPDinterval[2])
Rd.ps <- data.frame("Variable" = "Rd", "lambda"  = Rd.mod$lambdaMean, 
                        "lower HPD" = Rd.mod$lambdaHPDinterval[1], 
                        "upper HPD" = Rd.mod$lambdaHPDinterval[2])
log_gmax.ps <- data.frame("Variable" = "log_gmax", "lambda"  = log_gmax.mod$lambdaMean, 
                        "lower HPD" = log_gmax.mod$lambdaHPDinterval[1], 
                        "upper HPD" = log_gmax.mod$lambdaHPDinterval[2])

vars.ps <- rbind(PC1.ps, PC2.ps, PC3.ps, log_Ks.ps, log_Hv.ps, log_negMinWP_md.ps, log_negP50.ps, log_negPItlp.ps, log_Kl.ps, log_Hmax.ps, Rd.ps, log_gmax.ps)
vars.ps


########### HT ~ ENVIRONMENT UNIVARIATE #####################################################################################################

### log_Ks ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_PC1.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC1", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_PC1.mod$model)
log_Ks_PC1.mod$autocFix
log_Ks_PC1.mod$autocRan
log_Ks_PC1.mod$heidelFix
log_Ks_PC1.mod$heidelRan
log_Ks_PC1.mod$effSizeFix
log_Ks_PC1.mod$effSizeRan

# Lambda
plot(log_Ks_PC1.mod$lambda)
log_Ks_PC1.mod$lambdaMean
log_Ks_PC1.mod$lambdaHPDinterval

# Model results
summary(log_Ks_PC1.mod$model)

r2m <- R2mFun(log_Ks_PC1.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Ks_PC1.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Ks ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_PC2.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC2", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_PC2.mod$model)
log_Ks_PC2.mod$autocFix
log_Ks_PC2.mod$autocRan
log_Ks_PC2.mod$heidelFix
log_Ks_PC2.mod$heidelRan
log_Ks_PC2.mod$effSizeFix
log_Ks_PC2.mod$effSizeRan

# Lambda
plot(log_Ks_PC2.mod$lambda)
log_Ks_PC2.mod$lambdaMean
log_Ks_PC2.mod$lambdaHPDinterval

# Model results
summary(log_Ks_PC2.mod$model)

r2m <- R2mFun(log_Ks_PC2.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Ks_PC2.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Ks ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_PC3.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC3", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_PC3.mod$model)
log_Ks_PC3.mod$autocFix
log_Ks_PC3.mod$autocRan
log_Ks_PC3.mod$heidelFix
log_Ks_PC3.mod$heidelRan
log_Ks_PC3.mod$effSizeFix
log_Ks_PC3.mod$effSizeRan

# Lambda
plot(log_Ks_PC3.mod$lambda)
log_Ks_PC3.mod$lambdaMean
log_Ks_PC3.mod$lambdaHPDinterval

# Model results
summary(log_Ks_PC3.mod$model)

r2m <- R2mFun(log_Ks_PC3.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Ks_PC3.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negP50 ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_PC1.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC1", ran.frml = "~ animal + genus.Rand", 
                                  prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_PC1.mod$model)
log_negP50_PC1.mod$autocFix
log_negP50_PC1.mod$autocRan
log_negP50_PC1.mod$heidelFix
log_negP50_PC1.mod$heidelRan
log_negP50_PC1.mod$effSizeFix
log_negP50_PC1.mod$effSizeRan

# Lambda
plot(log_negP50_PC1.mod$lambda)
log_negP50_PC1.mod$lambdaMean
log_negP50_PC1.mod$lambdaHPDinterval

# Model results
summary(log_negP50_PC1.mod$model)

r2m <- R2mFun(log_negP50_PC1.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negP50_PC1.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negP50 ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_PC2.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC2", ran.frml = "~ animal + genus.Rand", 
                                  prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_PC2.mod$model)
log_negP50_PC2.mod$autocFix
log_negP50_PC2.mod$autocRan
log_negP50_PC2.mod$heidelFix
log_negP50_PC2.mod$heidelRan
log_negP50_PC2.mod$effSizeFix
log_negP50_PC2.mod$effSizeRan

# Lambda
plot(log_negP50_PC2.mod$lambda)
log_negP50_PC2.mod$lambdaMean
log_negP50_PC2.mod$lambdaHPDinterval

# Model results
summary(log_negP50_PC2.mod$model)


### log_negP50 ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_PC3.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC3", ran.frml = "~ animal + genus.Rand", 
                                  prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_PC3.mod$model)
log_negP50_PC3.mod$autocFix
log_negP50_PC3.mod$autocRan
log_negP50_PC3.mod$heidelFix
log_negP50_PC3.mod$heidelRan
log_negP50_PC3.mod$effSizeFix
log_negP50_PC3.mod$effSizeRan

# Lambda
plot(log_negP50_PC3.mod$lambda)
log_negP50_PC3.mod$lambdaMean
log_negP50_PC3.mod$lambdaHPDinterval

# Model results
summary(log_negP50_PC3.mod$model)

r2m <- R2mFun(log_negP50_PC3.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negP50_PC3.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negMinWP_md ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_PC1.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC1", ran.frml = "~ animal + genus.Rand", 
                                       prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_PC1.mod$model)
log_negMinWP_md_PC1.mod$autocFix
log_negMinWP_md_PC1.mod$autocRan
log_negMinWP_md_PC1.mod$heidelFix
log_negMinWP_md_PC1.mod$heidelRan
log_negMinWP_md_PC1.mod$effSizeFix
log_negMinWP_md_PC1.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_PC1.mod$lambda)
log_negMinWP_md_PC1.mod$lambdaMean
log_negMinWP_md_PC1.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_PC1.mod$model)

### log_negMinWP_md ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_PC2.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC2", ran.frml = "~ animal + genus.Rand", 
                                       prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_PC2.mod$model)
log_negMinWP_md_PC2.mod$autocFix
log_negMinWP_md_PC2.mod$autocRan
log_negMinWP_md_PC2.mod$heidelFix
log_negMinWP_md_PC2.mod$heidelRan
log_negMinWP_md_PC2.mod$effSizeFix
log_negMinWP_md_PC2.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_PC2.mod$lambda)
log_negMinWP_md_PC2.mod$lambdaMean
log_negMinWP_md_PC2.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_PC2.mod$model)

r2m <- R2mFun(log_negMinWP_md_PC2.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negMinWP_md_PC2.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negMinWP_md ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_PC3.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC3", ran.frml = "~ animal + genus.Rand", 
                                       prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_PC3.mod$model)
log_negMinWP_md_PC3.mod$autocFix
log_negMinWP_md_PC3.mod$autocRan
log_negMinWP_md_PC3.mod$heidelFix
log_negMinWP_md_PC3.mod$heidelRan
log_negMinWP_md_PC3.mod$effSizeFix
log_negMinWP_md_PC3.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_PC3.mod$lambda)
log_negMinWP_md_PC3.mod$lambdaMean
log_negMinWP_md_PC3.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_PC3.mod$model)

r2m <- R2mFun(log_negMinWP_md_PC3.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negMinWP_md_PC3.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Hv ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_PC1.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC1", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_PC1.mod$model)
log_Hv_PC1.mod$autocFix
log_Hv_PC1.mod$autocRan
log_Hv_PC1.mod$heidelFix
log_Hv_PC1.mod$heidelRan
log_Hv_PC1.mod$effSizeFix
log_Hv_PC1.mod$effSizeRan

# Lambda
plot(log_Hv_PC1.mod$lambda)
log_Hv_PC1.mod$lambdaMean
log_Hv_PC1.mod$lambdaHPDinterval

# Model results
summary(log_Hv_PC1.mod$model)

r2m <- R2mFun(log_Hv_PC1.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Hv_PC1.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Hv ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_PC2.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC2", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_PC2.mod$model)
log_Hv_PC2.mod$autocFix
log_Hv_PC2.mod$autocRan
log_Hv_PC2.mod$heidelFix
log_Hv_PC2.mod$heidelRan
log_Hv_PC2.mod$effSizeFix
log_Hv_PC2.mod$effSizeRan

# Lambda
plot(log_Hv_PC2.mod$lambda)
log_Hv_PC2.mod$lambdaMean
log_Hv_PC2.mod$lambdaHPDinterval

# Model results
summary(log_Hv_PC2.mod$model)

r2m <- R2mFun(log_Hv_PC2.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Hv_PC2.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Hv ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_PC3.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC3", ran.frml = "~ animal + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_PC3.mod$model)
log_Hv_PC3.mod$autocFix
log_Hv_PC3.mod$autocRan
log_Hv_PC3.mod$heidelFix
log_Hv_PC3.mod$heidelRan
log_Hv_PC3.mod$effSizeFix
log_Hv_PC3.mod$effSizeRan

# Lambda
plot(log_Hv_PC3.mod$lambda)
log_Hv_PC3.mod$lambdaMean
log_Hv_PC3.mod$lambdaHPDinterval

# Model results
summary(log_Hv_PC3.mod$model)

r2m <- R2mFun(log_Hv_PC3.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Hv_PC3.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

########### HT ~ HT ##############################################################################################################


### log_Ks ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_negP50.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ log_negP50", ran.frml = "~ animal + genus.Rand", 
                                     prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_log_negP50.mod$model)
log_Ks_log_negP50.mod$autocFix
log_Ks_log_negP50.mod$autocRan
log_Ks_log_negP50.mod$heidelFix
log_Ks_log_negP50.mod$heidelRan
log_Ks_log_negP50.mod$effSizeFix
log_Ks_log_negP50.mod$effSizeRan

# Lambda
plot(log_Ks_log_negP50.mod$lambda)
log_Ks_log_negP50.mod$lambdaMean
log_Ks_log_negP50.mod$lambdaHPDinterval

# Model results
summary(log_Ks_log_negP50.mod$model)

r2m <- R2mFun(log_Ks_log_negP50.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Ks_log_negP50.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negP50 ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_log_Ks.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ log_Ks", ran.frml = "~ animal + genus.Rand", 
                                     prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_log_Ks.mod$model)
log_negP50_log_Ks.mod$autocFix
log_negP50_log_Ks.mod$autocRan
log_negP50_log_Ks.mod$heidelFix
log_negP50_log_Ks.mod$heidelRan
log_negP50_log_Ks.mod$effSizeFix
log_negP50_log_Ks.mod$effSizeRan

# Lambda
plot(log_negP50_log_Ks.mod$lambda)
log_negP50_log_Ks.mod$lambdaMean
log_negP50_log_Ks.mod$lambdaHPDinterval

# Model results
summary(log_negP50_log_Ks.mod$model)

r2m <- R2mFun(log_negP50_log_Ks.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negP50_log_Ks.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negP50 ~ log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negP50_log_negMinWP_md.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ log_negMinWP_md", ran.frml = "~ animal + genus.Rand", 
                                              prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_log_negMinWP_md.mod$model)
log_negP50_log_negMinWP_md.mod$autocFix
log_negP50_log_negMinWP_md.mod$autocRan
log_negP50_log_negMinWP_md.mod$heidelFix
log_negP50_log_negMinWP_md.mod$heidelRan
log_negP50_log_negMinWP_md.mod$effSizeFix
log_negP50_log_negMinWP_md.mod$effSizeRan

# Lambda
plot(log_negP50_log_negMinWP_md.mod$lambda)
log_negP50_log_negMinWP_md.mod$lambdaMean
log_negP50_log_negMinWP_md.mod$lambdaHPDinterval

# Model results
summary(log_negP50_log_negMinWP_md.mod$model)


### log_negMinWP_md ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negMinWP_md_log_negP50.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ log_negP50", ran.frml = "~ animal + genus.Rand", 
                                              prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_log_negP50.mod$model)
log_negMinWP_md_log_negP50.mod$autocFix
log_negMinWP_md_log_negP50.mod$autocRan
log_negMinWP_md_log_negP50.mod$heidelFix
log_negMinWP_md_log_negP50.mod$heidelRan
log_negMinWP_md_log_negP50.mod$effSizeFix
log_negMinWP_md_log_negP50.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_log_negP50.mod$lambda)
log_negMinWP_md_log_negP50.mod$lambdaMean
log_negMinWP_md_log_negP50.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_log_negP50.mod$model)

r2m <- R2mFun(log_negMinWP_md_log_negP50.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negMinWP_md_log_negP50.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Hv ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_log_Ks.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ log_Ks", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_log_Ks.mod$model)
log_Hv_log_Ks.mod$autocFix
log_Hv_log_Ks.mod$autocRan
log_Hv_log_Ks.mod$heidelFix
log_Hv_log_Ks.mod$heidelRan
log_Hv_log_Ks.mod$effSizeFix
log_Hv_log_Ks.mod$effSizeRan

# Lambda
plot(log_Hv_log_Ks.mod$lambda)
log_Hv_log_Ks.mod$lambdaMean
log_Hv_log_Ks.mod$lambdaHPDinterval

# Model results
summary(log_Hv_log_Ks.mod$model)

r2m <- R2mFun(log_Hv_log_Ks.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Hv_log_Ks.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_Ks ~ log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_Hv.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ log_Hv", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_log_Hv.mod$model)
log_Ks_log_Hv.mod$autocFix
log_Ks_log_Hv.mod$autocRan
log_Ks_log_Hv.mod$heidelFix
log_Ks_log_Hv.mod$heidelRan
log_Ks_log_Hv.mod$effSizeFix
log_Ks_log_Hv.mod$effSizeRan

# Lambda
plot(log_Ks_log_Hv.mod$lambda)
log_Ks_log_Hv.mod$lambdaMean
log_Ks_log_Hv.mod$lambdaHPDinterval

# Model results
summary(log_Ks_log_Hv.mod$model)

r2m <- R2mFun(log_Ks_log_Hv.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Ks_log_Hv.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

############################# HT + HT ~ -1 ##############################################################################

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
meanKsP50.multiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1, 
                               random = ~us(trait):animal + us(trait):genus.Rand, 
                               rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                               prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsP50.multiMod)
plot(meanKsP50.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsP50.multiMod$Sol)
autocorr.diag(meanKsP50.multiMod$VCV)
# Convergence
heidel.diag(meanKsP50.multiMod$Sol)
heidel.diag(meanKsP50.multiMod$VCV)

## Phylogenetic signal

# log_Ks
physig_log_Ks <- meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_negP50
physig_log_negP50 <- meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)

## Inter-genus

# log_Ks
inter_log_Ks <- meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(inter_log_Ks)
HPDinterval(inter_log_Ks, 0.95)

# log_negP50
inter_log_negP50 <- meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(inter_log_negP50)
HPDinterval(inter_log_negP50, 0.95)

## Intra-genus

# log_Ks
intra_log_Ks <- meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(intra_log_Ks)
HPDinterval(intra_log_Ks, 0.95)

# log_negP50
intra_log_negP50 <- meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(intra_log_negP50)
HPDinterval(intra_log_negP50, 0.95)

## Evolutive correlation
cor_log_Ks_log_negP50.multiMod <- meanKsP50.multiMod$VCV[, "traitlog_negP50:traitlog_Ks.animal"]/
  sqrt(meanKsP50.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*meanKsP50.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(cor_log_Ks_log_negP50.multiMod)
HPDinterval(cor_log_Ks_log_negP50.multiMod, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
meanKsHv.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1, 
                              random = ~us(trait):animal + us(trait):genus.Rand, 
                              rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                              prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsHv.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsHv.multiMod$Sol)
autocorr.diag(meanKsHv.multiMod$VCV)
# Convergence
heidel.diag(meanKsHv.multiMod$Sol)
heidel.diag(meanKsHv.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks <- meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_Hv
physig_log_Hv <- meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Inter- genus

# log_Ks
intra_log_Ks <- meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(intra_log_Ks)
HPDinterval(intra_log_Ks, 0.95)

# log_Hv
intra_log_Hv <- meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"]/ 
  (meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(intra_log_Hv)
HPDinterval(intra_log_Hv, 0.95)

# Intra-genus

# log_Ks
intra_log_Ks <- meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsHv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(intra_log_Ks)
HPDinterval(intra_log_Ks, 0.95)

# log_Hv
intra_log_Hv <- meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"]/ 
  (meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(intra_log_Hv)
HPDinterval(intra_log_Hv, 0.95)

# Evolutive correlation
cor_log_Ks_log_Hv.multiMod <- meanKsHv.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(meanKsHv.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*meanKsHv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv.multiMod)
HPDinterval(cor_log_Ks_log_Hv.multiMod, 0.95)


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
meanMinWPP50.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1, 
                                  random = ~us(trait):animal + us(trait):genus.Rand, 
                                  rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                  prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanMinWPP50.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanMinWPP50.multiMod$Sol)
autocorr.diag(meanMinWPP50.multiMod$VCV)
# Convergence
heidel.diag(meanMinWPP50.multiMod$Sol)
heidel.diag(meanMinWPP50.multiMod$VCV)


## Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md <- meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_negP50
physig_log_negP50 <- meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)

# Inter-genus

# log_negMinWP_md
inter_log_negMinWP_md <- meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"]/ 
  (meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(inter_log_negMinWP_md)
HPDinterval(inter_log_negMinWP_md, 0.95)

# log_negP50
inter_log_negP50 <- meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(inter_log_negP50)
HPDinterval(inter_log_negP50, 0.95)


# Intra-genus

# log_negMinWP_md
intra_log_negMinWP_md <- meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]/ 
  (meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanMinWPP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(intra_log_negMinWP_md)
HPDinterval(intra_log_negMinWP_md, 0.95)

# log_negP50
intra_log_negP50 <- meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanMinWPP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(intra_log_negP50)
HPDinterval(intra_log_negP50, 0.95)

# Evolutive correlation
cor_log_negP50_log_negMinWP_md.multiMod <- meanMinWPP50.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negP50.animal"]/
  sqrt(meanMinWPP50.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*meanMinWPP50.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(corr.evo)
HPDinterval(cor_log_negP50_log_negMinWP_md.multiMod, 0.95)


### log_negMinWP_md + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
meanMinWPHv.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_Hv)~trait-1, 
                                  random = ~us(trait):animal + us(trait):genus.Rand, 
                                  rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                  prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanMinWPHv.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanMinWPHv.multiMod$Sol)
autocorr.diag(meanMinWPHv.multiMod$VCV)
# Convergence
heidel.diag(meanMinWPHv.multiMod$Sol)
heidel.diag(meanMinWPHv.multiMod$VCV)


## Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md <- meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_Hv
physig_log_Hv <- meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Inter-genus

# log_negMinWP_md
inter_log_negMinWP_md <- meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"]/ 
  (meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(inter_log_negMinWP_md)
HPDinterval(inter_log_negMinWP_md, 0.95)

# log_Hv
inter_log_Hv <- meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"]/ 
  (meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + meanMinWPHv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(inter_log_Hv)
HPDinterval(inter_log_Hv, 0.95)


# Intra-genus

# log_negMinWP_md
intra_log_negMinWP_md <- meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]/ 
  (meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanMinWPHv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(intra_log_negMinWP_md)
HPDinterval(intra_log_negMinWP_md, 0.95)

# log_Hv
intra_log_Hv <- meanMinWPP50.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"]/ 
  (meanMinWPP50.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanMinWPP50.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + meanMinWPP50.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(intra_log_Hv)
HPDinterval(intra_log_Hv, 0.95)

# Evolutive correlation
cor_log_Hv_log_negMinWP_md.multiMod <- meanMinWPHv.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_Hv.animal"]/
  sqrt(meanMinWPHv.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*meanMinWPHv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(corr.evo)
HPDinterval(cor_log_Hv_log_negMinWP_md.multiMod, 0.95)


### log_Hv + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Hv"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
meanHvP50.multiMod <- MCMCglmm(cbind(log_Hv, log_negP50)~trait-1, 
                                  random = ~us(trait):animal + us(trait):genus.Rand, 
                                  rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                  prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanHvP50.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanHvP50.multiMod$Sol)
autocorr.diag(meanHvP50.multiMod$VCV)
# Convergence
heidel.diag(meanHvP50.multiMod$Sol)
heidel.diag(meanHvP50.multiMod$VCV)


## Phylogenetic signal

# log_Hv
physig_log_Hv <- meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]/ 
  (meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"] + meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.genus.Rand"] + meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# log_negP50
physig_log_negP50 <- meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)

# Inter-genus

# log_Hv
inter_log_Hv <- meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.genus.Rand"]/ 
  (meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"] + meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.genus.Rand"] + meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.units"])
mean(inter_log_Hv)
HPDinterval(inter_log_Hv, 0.95)

# log_negP50
inter_log_negP50 <- meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(inter_log_negP50)
HPDinterval(inter_log_negP50, 0.95)


# Intra-genus

# log_Hv
intra_log_Hv <- meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.units"]/ 
  (meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"] + meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.genus.Rand"] + meanHvP50.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.units"])
mean(intra_log_Hv)
HPDinterval(intra_log_Hv, 0.95)

# log_negP50
intra_log_negP50 <- meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + meanHvP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(intra_log_negP50)
HPDinterval(intra_log_negP50, 0.95)

# Evolutive correlation
cor_log_negP50_log_Hv.multiMod <- meanHvP50.multiMod$VCV[, "traitlog_Hv:traitlog_negP50.animal"]/
  sqrt(meanHvP50.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]*meanHvP50.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(corr.evo)
HPDinterval(cor_log_negP50_log_Hv.multiMod, 0.95)


### log_Ks + log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negMinWP_md", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
meanKsMinWP.multiMod <- MCMCglmm(cbind(log_Ks, log_negMinWP_md)~trait-1, 
                               random = ~us(trait):animal + us(trait):genus.Rand, 
                               rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                               prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsMinWP.multiMod)
plot(meanKsMinWP.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsMinWP.multiMod$Sol)
autocorr.diag(meanKsMinWP.multiMod$VCV)
# Convergence
heidel.diag(meanKsMinWP.multiMod$Sol)
heidel.diag(meanKsMinWP.multiMod$VCV)

## Phylogenetic signal

# log_Ks
physig_log_Ks <- meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_negMinWP_md
physig_log_negMinWP_md <- meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

## Inter-genus

# log_Ks
inter_log_Ks <- meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(inter_log_Ks)
HPDinterval(inter_log_Ks, 0.95)

# log_negMinWP_md
inter_log_negMinWP_md <- meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"]/ 
  (meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(inter_log_negMinWP_md)
HPDinterval(inter_log_negMinWP_md, 0.95)

## Intra-genus

# log_Ks
intra_log_Ks <- meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + meanKsMinWP.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(intra_log_Ks)
HPDinterval(intra_log_Ks, 0.95)

# log_negMinWP_md
intra_log_negMinWP_md <- meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.units"]/ 
  (meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(intra_log_negMinWP_md)
HPDinterval(intra_log_negMinWP_md, 0.95)

## Evolutive correlation
cor_log_Ks_log_negMinWP_md.multiMod <- meanKsMinWP.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_Ks.animal"]/
  sqrt(meanKsMinWP.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*meanKsMinWP.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]) 
mean(cor_log_Ks_log_negMinWP_md.multiMod)
HPDinterval(cor_log_Ks_log_negMinWP_md.multiMod, 0.95)


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
names(data)[names(data) == "genus"] <- "animal"

# Model
log_Ks_log_negP50.multiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                       random = ~us(trait):animal + us(trait):genus.Rand, 
                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                       prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50.multiMod)
plot(log_Ks_log_negP50.multiMod)

r2m <- R2mFun(log_Ks_log_negP50.multiMod)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_Ks_log_negP50.multiMod)
r2c$R2c_mean
r2c$R2c_HPD

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50.multiMod$Sol)
autocorr.diag(log_Ks_log_negP50.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50.multiMod$Sol)
heidel.diag(log_Ks_log_negP50.multiMod$VCV)

## Phylogenetic signal

# log_Ks
physig_log_Ks <- log_Ks_log_negP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_negP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_negP50.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_negP50
physig_log_negP50 <- log_Ks_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (log_Ks_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_Ks_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)

## Evolutive correlation
cor_log_Ks_log_negP50_env.multiMod <- log_Ks_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_negP50.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_negP50.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(cor_log_Ks_log_negP50_env)
HPDinterval(cor_log_Ks_log_negP50_env, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_Ks_log_Hv.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                   random = ~us(trait):animal + us(trait):genus.Rand, 
                                   rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                   prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv.multiMod$Sol)
autocorr.diag(log_Ks_log_Hv.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv.multiMod$Sol)
heidel.diag(log_Ks_log_Hv.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks <- log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_Hv
physig_log_Hv <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Evolutive correlation
cor_log_Ks_log_Hv_env <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv_env)
HPDinterval(cor_log_Ks_log_Hv_env, 0.95)

### log_negMinWP_md + log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_negMinWP_md_log_Ks.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_Ks)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                                random = ~us(trait):animal + us(trait):genus.Rand, 
                                                rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                                prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_Ks.multiMod)
plot(log_negMinWP_md_log_Ks.multiMod)
## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_Ks.multiMod$Sol)
autocorr.diag(log_negMinWP_md_log_Ks.multiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_Ks.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_Ks.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md <- log_negMinWP_md_log_Ks.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_Ks.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_Ks.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_Ks
physig_log_Ks <- log_negMinWP_md_log_Ks.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_negMinWP_md_log_Ks.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"] + log_negMinWP_md_log_Ks.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# Evolutive correlation
cor_log_negMinWP_md_log_Ks_env <- log_negMinWP_md_log_Ks.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_Ks.animal"]/
  sqrt(log_negMinWP_md_log_Ks.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*log_negMinWP_md_log_Ks.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]) 
mean(cor_log_negMinWP_md_log_Ks_env)
HPDinterval(cor_log_negMinWP_md_log_Ks_env, 0.95)

### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_Ks_log_Hv.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                   random = ~us(trait):animal + us(trait):genus.Rand, 
                                   rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                   prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv.multiMod$Sol)
autocorr.diag(log_Ks_log_Hv.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv.multiMod$Sol)
heidel.diag(log_Ks_log_Hv.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks <- log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_Hv
physig_log_Hv <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Evolutive correlation
cor_log_Ks_log_Hv_env <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv_env)
HPDinterval(cor_log_Ks_log_Hv_env, 0.95)

### log_negMinWP_md + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_negMinWP_md_log_Hv.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_Hv)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                                random = ~us(trait):animal + us(trait):genus.Rand, 
                                                rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                                prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_Hv.multiMod)
plot(log_negMinWP_md_log_Hv.multiMod)
## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_Hv.multiMod$Sol)
autocorr.diag(log_negMinWP_md_log_Hv.multiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_Hv.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_Hv.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md <- log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_Hv
physig_log_Hv <- log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Evolutive correlation
cor_log_negMinWP_md_log_Hv_env <- log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_Hv.animal"]/
  sqrt(log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_negMinWP_md_log_Hv_env)
HPDinterval(cor_log_negMinWP_md_log_Hv_env, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_Ks_log_Hv.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                   random = ~us(trait):animal + us(trait):genus.Rand, 
                                   rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                   prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv.multiMod$Sol)
autocorr.diag(log_Ks_log_Hv.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv.multiMod$Sol)
heidel.diag(log_Ks_log_Hv.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks <- log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_Hv
physig_log_Hv <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Evolutive correlation
cor_log_Ks_log_Hv_env <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv_env)
HPDinterval(cor_log_Ks_log_Hv_env, 0.95)


### log_negMinWP_md + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_negMinWP_md_log_Hv.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_Hv)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, 
                                            random = ~us(trait):animal + us(trait):genus.Rand, 
                                            rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                            prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_Hv.multiMod)
plot(log_negMinWP_md_log_Hv.multiMod)
## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_Hv.multiMod$Sol)
autocorr.diag(log_negMinWP_md_log_Hv.multiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_Hv.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_Hv.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md <- log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_Hv
physig_log_Hv <- log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)

# Evolutive correlation
cor_log_negMinWP_md_log_Hv_env <- log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_Hv.animal"]/
  sqrt(log_negMinWP_md_log_Hv.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*log_negMinWP_md_log_Hv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_negMinWP_md_log_Hv_env)
HPDinterval(cor_log_negMinWP_md_log_Hv_env, 0.95)


#### SAVE IMAGE ####
save.image(file = "RData/Phylogenetic_linear_mixed_multi_models.RData")
