################################# PHYLOGENETIC LINEAR MIXED MODELS ##################################

# Pablo Sanchez
# 02-05-2019

#http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11

#install.packages("MCMCglmm")


### PACKAGES ###

library(ape)
library(geiger)
library(MCMCglmm)

### FUNCTIONS ####

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

mcmcglmmFun <- function(dta, fix.frml, ran.frml, prior, nitte = 50000, burn = 1000, thinn = 100){
  # Linear relationship between respnse and cofactor, with a random effect phylo corresponding the phylogenetic effect. The argument ginverse allows us to include a custom matrix for our random effect phylo, 
  mdl <- MCMCglmm(as.formula(fix.frml), random = as.formula(ran.frml),
                  family = "gaussian", ginverse = list(genus = inv.phylo$Ainv), prior = prior,
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
  lambda <- mdl$VCV[,'genus']/
    (mdl$VCV[,'genus'] + mdl$VCV[,'units'])
  rslts$lambda <- lambda
  rslts$lambdaMean <- mean(lambda)
  rslts$lambdaPost.mode <- posterior.mode(lambda)
  rslts$lambdaHPDinterval <- HPDinterval(lambda) # 95% credible interval 
  return(rslts)
}

### DATA ####
setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")
# RData
load("RData/Phylogenetic_linear_mixed_models.RData")

# Species datset
ht <- read.csv("Ht_Pablo/HydraTRY_Pablo_2019.csv", header = T)  # Includes three principal components

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

# Complete data
plmm.data <- completeFun(plmm.data, c("PC1", "Species"))

### PARAMETERS ####

## Matrix of phylogenetic correlation
# In order to include the phylogenetic effect we are going to calculate the inverse of the E matrix of phylogenetic correlation

inv.phylo <- inverseA(phylo, nodes="TIPS", scale=TRUE) # nodes = "ALL" includes ancestral nodes into the calculation

## Prior distributions for the fixed and random effects
# G for the random effects and R for the residual variance, corresponds to a inverse-Gamma distribution with shape and scape parameter equal to 0.01, which is relatively canonical. 
prior1 <- list(G = list(G1=list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002)) # we use the default prior for the fixed effects
prior2 <- list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002))) # two random effects, phylogeny and genus

priorMulti <- list(R=list(V=diag(2)/2,nu=2), G=list(G1=list(V=diag(2)/2,nu=2)))

########### HT ~ ENVIRONMENT #####################################################################################################
### log_Ks ~ environment ####

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Ks"))
length(data$Species)
# Model
log_Ks_env.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC1 + PC2 + PC3", ran.frml = "~ genus + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_env.mod$model)
log_Ks_env.mod$autocFix
log_Ks_env.mod$autocRan
log_Ks_env.mod$heidelFix
log_Ks_env.mod$heidelRan
log_Ks_env.mod$effSizeFix
log_Ks_env.mod$effSizeRan

# Lambda
plot(log_Ks_env.mod$lambda)
log_Ks_env.mod$lambdaMean
log_Ks_env.mod$lambdaHPDinterval

# Model results
summary(log_Ks_env.mod$model)

### log_negP50 ~ environment ####

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negP50"))
length(data$Species)
# Model
log_negP50_env.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC1 + PC2 + PC3", ran.frml = "~ genus + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_env.mod$model)
log_negP50_env.mod$autocFix
log_negP50_env.mod$autocRan
log_negP50_env.mod$heidelFix
log_negP50_env.mod$heidelRan
log_negP50_env.mod$effSizeFix
log_negP50_env.mod$effSizeRan

# Lambda
plot(log_negP50_env.mod$lambda)
log_negP50_env.mod$lambdaMean
log_negP50_env.mod$lambdaHPDinterval

# Model results
summary(log_negP50_env.mod$model)

### log_negMinWP_md ~ environment ####

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negMinWP_md"))
length(data$Species)
# Model
log_negMinWP_md_env.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC1 + PC2 + PC3", ran.frml = "~ genus + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_env.mod$model)
log_negMinWP_md_env.mod$autocFix
log_negMinWP_md_env.mod$autocRan
log_negMinWP_md_env.mod$heidelFix
log_negMinWP_md_env.mod$heidelRan
log_negMinWP_md_env.mod$effSizeFix
log_negMinWP_md_env.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_env.mod$lambda)
log_negMinWP_md_env.mod$lambdaMean
log_negMinWP_md_env.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_env.mod$model)

### log_Hv ~ environment ####

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Hv"))
length(data$Species)
# Model
log_Hv_env.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC1 + PC2 + PC3", ran.frml = "~ genus + genus.Rand", 
                              prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_env.mod$model)
log_Hv_env.mod$autocFix
log_Hv_env.mod$autocRan
log_Hv_env.mod$heidelFix
log_Hv_env.mod$heidelRan
log_Hv_env.mod$effSizeFix
log_Hv_env.mod$effSizeRan

# Lambda
plot(log_Hv_env.mod$lambda)
log_Hv_env.mod$lambdaMean
log_Hv_env.mod$lambdaHPDinterval

# Model results
summary(log_Hv_env.mod$model)

########### HT ~ HT ##############################################################################################################


### log_Ks ~ log_negP50 ####

# Complete data
Ks_P50.data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
length(Ks_P50.data$Species)

log_Ks_log_negP50.mod <- mcmcglmmFun(dta = Ks_P50.data, fix.frml = "log_Ks ~ log_negP50", ran.frml = "~ genus + genus.Rand", 
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

### log_negP50 ~ log_negMinWP_md ####

# Complete data
P50_MinWP.data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
length(P50_MinWP.data$Species)

log_negP50_log_negMinWP_md.mod <- mcmcglmmFun(dta = P50_MinWP.data, fix.frml = "log_negP50 ~ log_negMinWP_md", ran.frml = "~ genus + genus.Rand", 
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

### log_Hv ~ log_Ks ####

# Complete data
Hv_Ks.data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
length(Hv_Ks.data$Species)


log_Hv_log_Ks.mod <- mcmcglmmFun(dta = Hv_Ks.data, fix.frml = "log_Hv ~ log_Ks", ran.frml = "~ genus + genus.Rand", 
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


############################# HT + HT ##############################################################################

### log_Ks + log_negP50 ####

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_Ks_log_negP50.multiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, random = ~us(trait):animal, 
                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                       prior = priorMulti, data= data, nitt=260000,burnin=10000,thin=250)

summary(log_Ks_log_negP50.multiMod)

# Autocorrelation and convergence
heidel.diag(log_Ks_log_negP50.multiMod$Sol)
heidel.diag(log_Ks_log_negP50.multiMod$VCV)

# Phylogenetic signal

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


### log_Ks + log_Hv ####

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_Ks_log_Hv.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, random = ~us(trait):animal, 
                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                       prior = priorMulti, data= data, nitt=260000,burnin=10000,thin=250)

summary(log_Ks_log_Hv.multiMod)

# Autocorrelation and convergence
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
corr.evo <- log_Ks_log_Hv.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(corr.evo)
HPDinterval(corr.evo, 0.95)


### log_negMinWP_md + log_negP50 ####

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
names(data)[names(data) == "genus"] <- "animal"

# Model
log_negMinWP_md_log_negP50.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1 + trait:PC1 + trait:PC2 + trait:PC3, random = ~us(trait):animal, 
                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                       prior = priorMulti, data= data, nitt=260000,burnin=10000,thin=250)

summary(log_negMinWP_md_log_negP50.multiMod)

# Autocorrelation and convergence
heidel.diag(log_negMinWP_md_log_negP50.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_negP50.multiMod$VCV)

# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md <- log_negMinWP_md_log_negP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_negP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_negP50.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_negP50
physig_log_negP50 <- log_negMinWP_md_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (log_negMinWP_md_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_negMinWP_md_log_negP50.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)

#### SAVE IMAGE ####

#save.image(file = "RData/Phylogenetic_linear_mixed_models.RData")
