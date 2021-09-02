################################# COVARIATION STRUCTURES ##################################

# Pablo Sanchez
# 02-05-2019

#http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11

### PACKAGES ###
library(ape)
library(MCMCglmm)

setwd("~/Pablo/data")

# RDATA
load("RData/not_used/cv_structure.RData")

load("RData/Phylogenetic_linear_mixed_multi_models_groups.RData")
load("RData/Phylogenetic_linear_mixed_multi_models.RData")

# Functions
multiModelsIndexFun <- function(mdl, v1, v2) {
  if(length(v2) == 1){
    vars <- c("PC1", "PC2", "PC3", v1, v2)
  } else {
    vars <- c("PC1", "PC2", "PC3", v1)
  }
  rslts <- list()
  n <- 0
  e <- 0
  phylo.df <- data.frame()
  corr.df <- data.frame()
  for(var1 in vars){
  n <- n + 1
  phylo.df[n, "variable"] <- var1
  ## Phylogenetic signal
  lambda <- mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")]/
    (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")] + 
       mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
  phylo.df[n, "lambdaMean"] <- mean(lambda)
  phylo.df[n, "lambda_lowerHPD"] <- HPDinterval(lambda)[1] # 95% credible interval 
  phylo.df[n, "lambda_upperHPD"] <- HPDinterval(lambda)[2]
  ## inter-genus
  inter <- mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")]/
    (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")] + 
       mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
  phylo.df[n, "interMean"] <- mean(inter)
  phylo.df[n, "inter_lowerHDP"] <- HPDinterval(inter)[1] # 95% credible interval 
  phylo.df[n, "inter_upperHDP"] <- HPDinterval(inter)[2]
  ## Intra-genus
  intra <- mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")]/
    (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")] + 
       mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
  phylo.df[n, "intraMean"] <- mean(intra)
  phylo.df[n, "intra_lowerHDP"] <- HPDinterval(intra)[1] # 95% credible interval 
  phylo.df[n, "intra_lowerHDP"] <- HPDinterval(intra)[2]
  rslts$phylo_var <-   phylo.df
  # Correlated evolution
 for(var2 in vars){
 e <- e + 1
 corr.df[e, "var1"] <- var1
 corr.df[e, "var2"] <- var2
 evoCor  <- mdl$VCV[, paste0("trait", var1, ":trait", var2, ".animal")]/
    sqrt(mdl$VCV[, paste0("trait", var2, ":trait", var2, ".animal")] * mdl$VCV[, paste0("trait", var1, ":trait", var1, ".animal")]) 
 corr.df[e, "meanCor"]  <- mean(evoCor)
 corr.df[e, "cor_lowerHDP"]  <- HPDinterval(evoCor, 0.95)[1]
 corr.df[e, "cor_upperHDP"]  <- HPDinterval(evoCor, 0.95)[2]
 rslts$evoCor <- corr.df
  }
  }
  return(rslts)
}

# Model parameters

priorMultiCVstr <- list(
  R=list(V=diag(5)*(0.002/1.002),nu=1.002), 
  G=list(G1=list(V=diag(5)*(0.002/1.002),nu=1.002),
         G2 = list(V=diag(5)/2,nu=2)))

priorMultiPCsCVstr <- list(
  R=list(V=diag(4)*(0.002/1.002),nu=1.002), 
  G=list(G1=list(V=diag(4)*(0.002/1.002),nu=1.002),
         G2 = list(V=diag(4)/2,nu=2)))

nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

#### WITHOUT GROUP AS FIXED EFFECT ###########################################################################################

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks", "PC1"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50_env_cvStr<- MCMCglmm(cbind(log_Ks, log_negP50, PC1, PC2, PC3)~trait-1, 
                                          random = ~us(trait):animal + us(trait):genus.Rand, 
                                          rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                          prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50_env_cvStr)
plot(log_Ks_log_negP50_env_cvStr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50_env_cvStr$Sol)
autocorr.diag(log_Ks_log_negP50_env_cvStr$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50_env_cvStr$Sol)
heidel.diag(log_Ks_log_negP50_env_cvStr$VCV)

## Results

log_Ks_log_negP50_env_Cstr.res <- multiModelsIndexFun(log_Ks_log_negP50_env_cvStr, v1 = "log_Ks", v2 = "log_negP50")

log_Ks_log_negP50_env_Cstr.res$phylo_var
log_Ks_log_negP50_env_Cstr.res$evoCor


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_Hv_env_cvStr <- MCMCglmm(cbind(log_Ks, log_Hv, PC1, PC2, PC3)~trait-1, 
                                    random = ~us(trait):animal + us(trait):genus.Rand, 
                                    rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                    prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv_env_cvStr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv_env_cvStr$Sol)
autocorr.diag(log_Ks_log_Hv_env_cvStr$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv_env_cvStr$Sol)
heidel.diag(log_Ks_log_Hv_env_cvStr$VCV)


## Results
log_Hv_log_Ks_env_Cstr.res <- multiModelsIndexFun(log_Ks_log_Hv_env_cvStr, v1 = "log_Hv", v2 = "log_Ks")
log_Hv_log_Ks_env_Cstr.res$phylo_var
log_Hv_log_Ks_env_Cstr.res$evoCor


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_negP50_env_CVstr <- MCMCglmm(cbind(log_negMinWP_md, log_negP50, PC1, PC2, PC3)~trait-1, 
                                             random = ~us(trait):animal + us(trait):genus.Rand, 
                                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                             prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_negP50_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_negP50_env_CVstr$Sol)
autocorr.diag(log_negMinWP_md_log_negP50_env_CVstr$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_negP50_env_CVstr$Sol)
heidel.diag(log_negMinWP_md_log_negP50_env_CVstr$VCV)


## Results
log_negMinWP_md_log_negP50_env_Cstr.res <- multiModelsIndexFun(log_negMinWP_md_log_negP50_env_CVstr, v1 = "log_negMinWP_md", v2 = "log_negP50")

log_negMinWP_md_log_negP50_env_Cstr.res$phylo_var
log_negMinWP_md_log_negP50_env_Cstr.res$evoCor


### log_Hv + log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negMinWP_md", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_Hv_env_CVstr <- MCMCglmm(cbind(log_negMinWP_md, log_Hv, PC1, PC2, PC3)~trait-1, 
                                         random = ~us(trait):animal + us(trait):genus.Rand, 
                                         rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                         prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_Hv_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_Hv_env_CVstr$Sol)
autocorr.diag(log_negMinWP_md_log_Hv_env_CVstr$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_Hv_env_CVstr$Sol)
heidel.diag(log_negMinWP_md_log_Hv_env_CVstr$VCV)

## Results
log_negMinWP_md_log_Hv_env_Cstr.res <- multiModelsIndexFun(log_negMinWP_md_log_Hv_env_CVstr, v1 = "log_negMinWP_md", v2 = "log_Hv")

log_negMinWP_md_log_Hv_env_Cstr.res$phylo_var
log_negMinWP_md_log_Hv_env_Cstr.res$evoCor


### log_Hv + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negP50", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negP50_log_Hv_env_CVstr <- MCMCglmm(cbind(log_negP50, log_Hv, PC1, PC2, PC3)~trait-1, 
                                    random = ~us(trait):animal + us(trait):genus.Rand, 
                                    rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                    prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negP50_log_Hv_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negP50_log_Hv_env_CVstr$Sol)
autocorr.diag(log_negP50_log_Hv_env_CVstr$VCV)
# Convergence
heidel.diag(log_negP50_log_Hv_env_CVstr$Sol)
heidel.diag(log_negP50_log_Hv_env_CVstr$VCV)

## Results
log_negP50_log_Hv_env_Cstr.res <- multiModelsIndexFun(log_negP50_log_Hv_env_CVstr, v1 = "log_negP50", v2 = "log_Hv")

log_negP50_log_Hv_env_Cstr.res$phylo_var
log_negP50_log_Hv_env_Cstr.res$evoCor


#### WITH GROUP AS FIXED EFFECT #############################################################################################

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks", "PC1"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50_env_gr_cvStr<- MCMCglmm(cbind(log_Ks, log_negP50, PC1, PC2, PC3)~trait-1 + trait:group, 
                                                random = ~us(trait):animal + us(trait):genus.Rand, 
                                                rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                                prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50_env_gr_cvStr)
plot(log_Ks_log_negP50_env_gr_cvStr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50_env_gr_cvStr$Sol)
autocorr.diag(log_Ks_log_negP50_env_gr_cvStr$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50_env_gr_cvStr$Sol)
heidel.diag(log_Ks_log_negP50_env_gr_cvStr$VCV)

## Results
log_Ks_log_negP50_env_gr_Cstr.res <- multiModelsIndexFun(log_Ks_log_negP50_env_gr_cvStr, v1 = "log_Ks", v2 = "log_negP50")
log_Ks_log_negP50_env_gr_Cstr.res$phylo_var
log_Ks_log_negP50_env_gr_Cstr.res$evoCor


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_Hv_env_gr_cvStr <- MCMCglmm(cbind(log_Ks, log_Hv, PC1, PC2, PC3)~trait-1 + trait:group, 
                                    random = ~us(trait):animal + us(trait):genus.Rand, 
                                    rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                    prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv_env_gr_cvStr)
plot(log_Ks_log_Hv_env_gr_cvStr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv_env_gr_cvStr$Sol)
autocorr.diag(log_Ks_log_Hv_env_gr_cvStr$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv_env_gr_cvStr$Sol)
heidel.diag(log_Ks_log_Hv_env_gr_cvStr$VCV)

## Results
log_Hv_log_Ks_env_gr_Cstr.res <- multiModelsIndexFun(log_Ks_log_Hv_env_gr_cvStr, v1 = "log_Hv", v2 = "log_Ks")
log_Hv_log_Ks_env_gr_Cstr.res$phylo_var
log_Hv_log_Ks_env_gr_Cstr.res$evoCor


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_negP50_env_gr_CVstr <- MCMCglmm(cbind(log_negMinWP_md, log_negP50, PC1, PC2, PC3)~trait-1 + trait:group, 
                                             random = ~us(trait):animal + us(trait):genus.Rand, 
                                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                             prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_negP50_env_gr_CVstr)
plot(log_negMinWP_md_log_negP50_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_negP50_env_gr_CVstr$Sol)
autocorr.diag(log_negMinWP_md_log_negP50_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_negP50_env_gr_CVstr$Sol)
heidel.diag(log_negMinWP_md_log_negP50_env_gr_CVstr$VCV)


## Results
log_negMinWP_md_log_negP50_env_gr_Cstr.res <- multiModelsIndexFun(log_negMinWP_md_log_negP50_env_gr_CVstr, v1 = "log_negMinWP_md", v2 = "log_negP50")
log_negMinWP_md_log_negP50_env_gr_Cstr.res$phylo_var
log_negMinWP_md_log_negP50_env_gr_Cstr.res$evoCor


### log_Hv + log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negMinWP_md", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_Hv_env_gr_CVstr <- MCMCglmm(cbind(log_negMinWP_md, log_Hv, PC1, PC2, PC3)~trait-1 + trait:group, 
                                             random = ~us(trait):animal + us(trait):genus.Rand, 
                                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                             prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_Hv_env_gr_CVstr)
plot(log_negMinWP_md_log_Hv_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_Hv_env_gr_CVstr$Sol)
autocorr.diag(log_negMinWP_md_log_Hv_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_Hv_env_gr_CVstr$Sol)
heidel.diag(log_negMinWP_md_log_Hv_env_gr_CVstr$VCV)

## Results
log_negMinWP_md_log_Hv_env_gr_Cstr.res <- multiModelsIndexFun(log_negMinWP_md_log_Hv_env_gr_CVstr, v1 = "log_negMinWP_md", v2 = "log_Hv")
log_negMinWP_md_log_Hv_env_gr_Cstr.res$phylo_var
log_negMinWP_md_log_Hv_env_gr_Cstr.res$evoCor


### log_Hv + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negP50", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negP50_log_Hv_env_gr_CVstr <- MCMCglmm(cbind(log_negP50, log_Hv, PC1, PC2, PC3)~trait-1 + trait:group, 
                                         random = ~us(trait):animal + us(trait):genus.Rand, 
                                         rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                         prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negP50_log_Hv_env_gr_CVstr)
plot(log_negP50_log_Hv_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negP50_log_Hv_env_gr_CVstr$Sol)
autocorr.diag(log_negP50_log_Hv_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_negP50_log_Hv_env_gr_CVstr$Sol)
heidel.diag(log_negP50_log_Hv_env_gr_CVstr$VCV)

## Results
log_negP50_log_Hv_env_gr_Cstr.res <- multiModelsIndexFun(log_negP50_log_Hv_env_gr_CVstr, v1 = "log_negP50", v2 = "log_Hv")
log_negP50_log_Hv_env_gr_Cstr.res$phylo_var
log_negP50_log_Hv_env_gr_Cstr.res$evoCor


#### SAME DATASET #######################################################################################################

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks", "PC1", "log_Hv", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50_env_cvStr<- MCMCglmm(cbind(log_Ks, log_negP50, PC1, PC2, PC3)~trait-1, 
                                       random = ~us(trait):animal + us(trait):genus.Rand, 
                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                       prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50_env_cvStr)
plot(log_Ks_log_negP50_env_cvStr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50_env_cvStr$Sol)
autocorr.diag(log_Ks_log_negP50_env_cvStr$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50_env_cvStr$Sol)
heidel.diag(log_Ks_log_negP50_env_cvStr$VCV)

## Results

log_Ks_log_negP50_env_Cstr.res <- multiModelsIndexFun(log_Ks_log_negP50_env_cvStr, v1 = "log_Ks", v2 = "log_negP50")

log_Ks_log_negP50_env_Cstr.res$phylo_var
log_Ks_log_negP50_env_Cstr.res$evoCor


### log_Ks + log_Hv ###

# Model
log_Ks_log_Hv_env_cvStr <- MCMCglmm(cbind(log_Ks, log_Hv, PC1, PC2, PC3)~trait-1, 
                                    random = ~us(trait):animal + us(trait):genus.Rand, 
                                    rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian", "gaussian"),
                                    prior = priorMultiCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv_env_cvStr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv_env_cvStr$Sol)
autocorr.diag(log_Ks_log_Hv_env_cvStr$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv_env_cvStr$Sol)
heidel.diag(log_Ks_log_Hv_env_cvStr$VCV)


## Results
log_Hv_log_Ks_env_Cstr.res <- multiModelsIndexFun(log_Ks_log_Hv_env_cvStr, v1 = "log_Hv", v2 = "log_Ks")
log_Hv_log_Ks_env_Cstr.res$phylo_var
log_Hv_log_Ks_env_Cstr.res$evoCor

### ONLY PCS ###########################################################################################################

#### WHITOUT GROUP #####################################################################################################

### log_Ks + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_env_CVstr <- MCMCglmm(cbind(log_Ks, PC1, PC2, PC3)~trait-1, 
                                        random = ~us(trait):animal + us(trait):genus.Rand, 
                                        rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                                        prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_env_CVstr$Sol)
autocorr.diag(log_Ks_env_CVstr$VCV)
# Convergence
heidel.diag(log_Ks_env_CVstr$Sol)
heidel.diag(log_Ks_env_CVstr$VCV)

## Results
log_Ks_env_Cstr.res <- multiModelsIndexFun(log_Ks_env_CVstr, v1 = "log_Ks", v2 = NULL)

log_Ks_env_Cstr.res$phylo_var
log_Ks_env_Cstr.res$evoCor


### log_negP50 + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negP50_env_CVstr <- MCMCglmm(cbind(log_negP50, PC1, PC2, PC3)~trait-1, 
                             random = ~us(trait):animal + us(trait):genus.Rand, 
                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                             prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negP50_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negP50_env_CVstr$Sol)
autocorr.diag(log_negP50_env_CVstr$VCV)
# Convergence
heidel.diag(log_negP50_env_CVstr$Sol)
heidel.diag(log_negP50_env_CVstr$VCV)

## Results
log_negP50_env_Cstr.res <- multiModelsIndexFun(log_negP50_env_CVstr, v1 = "log_negP50", v2 = NULL)

log_negP50_env_Cstr.res$phylo_var
log_negP50_env_Cstr.res$evoCor

### log_negMinWP_md + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_negMinWP_md", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_env_CVstr <- MCMCglmm(cbind(log_negMinWP_md, PC1, PC2, PC3)~trait-1, 
                             random = ~us(trait):animal + us(trait):genus.Rand, 
                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                             prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_env_CVstr$Sol)
autocorr.diag(log_negMinWP_md_env_CVstr$VCV)
# Convergence
heidel.diag(log_negMinWP_md_env_CVstr$Sol)
heidel.diag(log_negMinWP_md_env_CVstr$VCV)

## Results
log_negMinWP_md_env_Cstr.res <- multiModelsIndexFun(log_negMinWP_md_env_CVstr, v1 = "log_negMinWP_md", v2 = NULL)

log_negMinWP_md_env_Cstr.res$phylo_var
log_negMinWP_md_env_Cstr.res$evoCor

### log_Hv + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Hv_env_CVstr <- MCMCglmm(cbind(log_Hv, PC1, PC2, PC3)~trait-1, 
                             random = ~us(trait):animal + us(trait):genus.Rand, 
                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                             prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Hv_env_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Hv_env_CVstr$Sol)
autocorr.diag(log_Hv_env_CVstr$VCV)
# Convergence
heidel.diag(log_Hv_env_CVstr$Sol)
heidel.diag(log_Hv_env_CVstr$VCV)

## Results
log_Hv_env_Cstr.res <- multiModelsIndexFun(log_Hv_env_CVstr, v1 = "log_Hv", v2 = NULL)

log_Hv_env_Cstr.res$phylo_var
log_Hv_env_Cstr.res$evoCor


#### WHITh GROUP #####################################################################################################

### log_Ks + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_env_gr_CVstr <- MCMCglmm(cbind(log_Ks, PC1, PC2, PC3)~trait-1 + trait:group, 
                             random = ~us(trait):animal + us(trait):genus.Rand, 
                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                             prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_env_gr_CVstr$Sol)
autocorr.diag(log_Ks_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_Ks_env_gr_CVstr$Sol)
heidel.diag(log_Ks_env_gr_CVstr$VCV)

## Results
log_Ks_env_gr_Cstr.res <- multiModelsIndexFun(log_Ks_env_gr_CVstr, v1 = "log_Ks", v2 = NULL)

log_Ks_env_gr_Cstr.res$phylo_var
log_Ks_env_gr_Cstr.res$evoCor


### log_negP50 + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negP50_env_gr_CVstr <- MCMCglmm(cbind(log_negP50, PC1, PC2, PC3)~trait-1 + trait:group, 
                                 random = ~us(trait):animal + us(trait):genus.Rand, 
                                 rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                                 prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negP50_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negP50_env_gr_CVstr$Sol)
autocorr.diag(log_negP50_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_negP50_env_gr_CVstr$Sol)
heidel.diag(log_negP50_env_gr_CVstr$VCV)

## Results
log_negP50_env_gr_Cstr.res <- multiModelsIndexFun(log_negP50_env_gr_CVstr, v1 = "log_negP50", v2 = NULL)

log_negP50_env_gr_Cstr.res$phylo_var
log_negP50_env_gr_Cstr.res$evoCor

### log_negMinWP_md + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_negMinWP_md", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_env_gr_CVstr <- MCMCglmm(cbind(log_negMinWP_md, PC1, PC2, PC3)~trait-1 + trait:group, 
                                      random = ~us(trait):animal + us(trait):genus.Rand, 
                                      rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                                      prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_env_gr_CVstr$Sol)
autocorr.diag(log_negMinWP_md_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_negMinWP_md_env_gr_CVstr$Sol)
heidel.diag(log_negMinWP_md_env_gr_CVstr$VCV)

## Results
log_negMinWP_md_env_gr_Cstr.res <- multiModelsIndexFun(log_negMinWP_md_env_gr_CVstr, v1 = "log_negMinWP_md", v2 = NULL)

log_negMinWP_md_env_gr_Cstr.res$phylo_var
log_negMinWP_md_env_gr_Cstr.res$evoCor

### log_Hv + PCs ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "PC1"))

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Hv_env_gr_CVstr <- MCMCglmm(cbind(log_Hv, PC1, PC2, PC3)~trait-1 + trait:group, 
                             random = ~us(trait):animal + us(trait):genus.Rand, 
                             rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian", "gaussian", "gaussian"),
                             prior = priorMultiPCsCVstr, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Hv_env_gr_CVstr)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Hv_env_gr_CVstr$Sol)
autocorr.diag(log_Hv_env_gr_CVstr$VCV)
# Convergence
heidel.diag(log_Hv_env_gr_CVstr$Sol)
heidel.diag(log_Hv_env_gr_CVstr$VCV)

## Results
log_Hv_env_gr_Cstr.res <- multiModelsIndexFun(log_Hv_env_gr_CVstr, v1 = "log_Hv", v2 = NULL)

log_Hv_env_gr_Cstr.res$phylo_var
log_Hv_env_gr_Cstr.res$evoCor


#### SAVE IMAGE ####
save.image(file = "RData/cv_structure.RData")
