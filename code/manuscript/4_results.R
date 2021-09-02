#### MODELS RESULTS ###################################################################################################################

# P Sanchez-Martinez

print("Running results...")
remove(list = ls())
source("code/manuscript/1_init.R")

print("loading multivariate models results...")
load("data/RData/plmm_multivariate_NP_results_diagn.RData")
load("data/RData/plmm_multivariate_results_diagn.RData")

print("loading univariate models results...")
load("data/RData/plmm_univariate_NP_results_diagn.RData")
load("data/RData/plmm_univariate_results_diagn.RData")

### FUNCTIONS ###

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
  rslts$R2m_upperHPD <- HPDinterval(R2m)[1]
  rslts$R2m_lowerHPD <- HPDinterval(R2m)[2]
  rslts <- as.data.frame(rslts)
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
  rslts$R2c_upperHPD <- HPDinterval(R2c)[1]
  rslts$R2c_lowerHPD <- HPDinterval(R2c)[2]
  rslts <- as.data.frame(rslts)
  return(rslts)
}

# No-phylogenetic R2GLMM(m) - marginal R2GLMM (eq. 26, 29 and 30 Nakawgaway 2013)
R2mNPFun <- function(mdl){
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
  rslts$R2m_upperHPD <- HPDinterval(R2m)[1]
  rslts$R2m_lowerHPD <- HPDinterval(R2m)[2]
  return(rslts)
}

# No-phylogenetic R2GLMM(c) - conditional R2GLMM for full model (equ. 30 Nakagawa 2013)
R2cNPFun <- function(mdl){
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
  rslts$R2c_upperHPD <- HPDinterval(R2c)[1]
  rslts$R2c_lowerHPD <- HPDinterval(R2c)[2]
  rslts <- as.data.frame(rslts)
  return(rslts)
}

# R2 results function
uniModR2Fun <- function(df.res){
  r2.df <- data.frame(model = rep(NA, length(df.res)))
  r2.df$R2m_mean <- numeric(length(df.res))
  r2.df$R2m_lowerHPD <- numeric(length(df.res))
  r2.df$R2m_upperHPD <- numeric(length(df.res))
  r2.df$R2c_mean <- numeric(length(df.res))
  r2.df$R2c_lowerHPD <- numeric(length(df.res))
  r2.df$R2c_upperHPD <- numeric(length(df.res))
  for(i in 1:length(df.res)){
    mod.name <- df.res[[i]]$"name"
    if(is.na(str_extract(mod.name, "NP"))){
    mod <- df.res[[i]]$"model"
    r2.df$model[i]  <- mod.name
    r2.df$R2m_mean[i] <- as.numeric(unlist(R2mFun(mod)[[2]]))
    r2.df$R2m_lowerHPD[i] <- as.numeric(unlist(R2mFun(mod)[[3]]))
    r2.df$R2m_upperHPD[i] <- as.numeric(unlist(R2mFun(mod)[[4]]))
    r2.df$R2c_mean[i] <-  as.numeric(unlist(R2cFun(mod)[[2]]))
    r2.df$R2c_lowerHPD[i] <- as.numeric(unlist(R2cFun(mod)[[3]]))
    r2.df$R2c_upperHPD[i] <- as.numeric(unlist(R2cFun(mod)[[4]]))
    }
    if(!is.na(str_extract(mod.name, "NP"))){
      mod <- df.res[[i]]$"model"
      r2.df$model[i]  <- mod.name
      r2.df$R2m_mean[i] <- as.numeric(unlist(R2mNPFun(mod)[[2]]))
      r2.df$R2m_lowerHPD[i] <- as.numeric(unlist(R2mNPFun(mod)[[3]]))
      r2.df$R2m_upperHPD[i] <- as.numeric(unlist(R2mNPFun(mod)[[4]]))
      r2.df$R2c_mean[i] <-  as.numeric(unlist(R2cNPFun(mod)[[2]]))
      r2.df$R2c_lowerHPD[i] <- as.numeric(unlist(R2cNPFun(mod)[[3]]))
      r2.df$R2c_upperHPD[i] <- as.numeric(unlist(R2cNPFun(mod)[[4]]))
    }
  }
  r2.df <- cbind("model" = r2.df$model, round(r2.df[, c(2:7)], 3))
  return(r2.df)
}

# DIC
DICFun <- function(mdls, mod.res){
  mdls$DIC <- NA
  for(i in 1:length(mdls$type)){
    mdls$DIC[i] <- mod.res[[i]]$DIC
  }
  mdls <- arrange(mdls,variables, desc(DIC))
  return(mdls)
}

# Best models
bestFun <- function(mdls, mdlsType = "Ht_Env_1"){
  bestMdls <- data.frame()
  mdls <- mdls[which(mdls$type == mdlsType), ]
  vars <- unique(mdls$variables)
  for(var in vars){
  mdlGroup <- mdls[which(mdls$variables == var), ]
  DICmin <- min(mdlGroup$DIC)
  bestMdl <- mdlGroup[which(mdlGroup$DIC < (DICmin + 4)), ]
  if(length(bestMdl$type) > 1){
    if(any(!is.na(str_extract(as.character(bestMdl$group), "NG")))){
      bestMdl <- bestMdl[which(bestMdl$group == "NG"), ]
      bestMdls <- rbind(bestMdls, bestMdl)
    } else{
      if(any(!is.na(str_extract(as.character(bestMdl$group), "GI")))){
        bestMdl <- bestMdl[which(bestMdl$group == "G"), ]
        bestMdls <- rbind(bestMdls, bestMdl)
      }
    }
  } else{
    bestMdls <- rbind(bestMdls, bestMdl) 
  }
  }
  return(bestMdls)
}

# Phylogenetic indexes

# Uni-models phylogenetic signal
uniModelsIndexFun <- function() {
  rslts <- list()
  phylo.df <- data.frame()
  for(n in 1:length(uniMdls$type)){ 
    var1 <- as.vector(unlist(strsplit(as.character(uniMdls[n, "variables"]), ", ")))[1]
    mdlName <- uniMdls[n, "name"]
    mdl <- uniMod.res[[mdlName]]$model
    nObs <- uniMod.res[[mdlName]]$NumObs
    if(is.na(str_extract(mdlName, "NP"))){
        phylo.df[n, "model"] <- mdlName
        phylo.df[n, "variable"] <- var1
        phylo.df[n, "N"] <- nObs
        ## Phylogenetic signal
        lambda <- mdl$VCV[,"animal"]/ (mdl$VCV[, "animal"] + mdl$VCV[, "genus.Rand"] + mdl$VCV[, "units"])
        pd <- p_direction(lambda)
        pd <- as.numeric(pd$pd)
        pval <- 2*(1 - pd)
        phylo.df[n, "lambdaMean"] <- mean(lambda)
        phylo.df[n, "lambda_lowerHPD"] <- HPDinterval(lambda)[1] # 95% credible interval 
        phylo.df[n, "lambda_upperHPD"] <- HPDinterval(lambda)[2]
        phylo.df[n, "lambda_p_value"] <- pval
        ## inter-genus
        inter <- mdl$VCV[, "genus.Rand"]/(mdl$VCV[, "animal"] + mdl$VCV[, "genus.Rand"] + mdl$VCV[, "units"])
        pd <- p_direction(inter)
        pd <- as.numeric(pd$pd)
        pval <- 2*(1 - pd)
        phylo.df[n, "interMean"] <- mean(inter)
        phylo.df[n, "inter_lowerHDP"] <- HPDinterval(inter)[1] # 95% credible interval 
        phylo.df[n, "inter_upperHDP"] <- HPDinterval(inter)[2]
        phylo.df[n, "inter_p_value"] <- pval
        ## Intra-genus
        intra <- mdl$VCV[, "units"]/(mdl$VCV[, "animal"] + mdl$VCV[, "genus.Rand"] + mdl$VCV[, "units"])
        phylo.df[n, "intraMean"] <- mean(intra)
        pd <- p_direction(intra)
        pd <- as.numeric(pd$pd)
        pval <- 2*(1 - pd)
        phylo.df[n, "intra_lowerHDP"] <- HPDinterval(intra)[1] # 95% credible interval 
        phylo.df[n, "intra_upperHDP"] <- HPDinterval(intra)[2]
        phylo.df[n, "intra_p_value"] <- pval
    }
  }
  # Order results
  phylo.df <- arrange(phylo.df, model)
  phylo.df <- unique(phylo.df)
  rslts$phylo_var <-   phylo.df
  return(rslts)
}

# Non-phylogenetic Uni-models phylogenetic signal
NP_uniModelsIndexFun <- function() {
  rslts <- list()
  phylo.df <- data.frame()
  for(n in 1:length(uniMdls$type)){ 
    vars <- as.vector(unlist(strsplit(as.character(uniMdls[n, "variables"]), ", ")))[1]
    mdl <- NP_uniMod.res[[n]]$"model"
    mdlName <- NP_uniMod.res[[n]]$"name"
    if(!is.na(str_extract(mdlName, "NP"))){
      for(var1 in vars){
        phylo.df[n, "model"] <- mdlName
        phylo.df[n, "variable"] <- var1
        ## inter-genus
        inter <- mdl$VCV[, "genus.Rand"]/(mdl$VCV[, "genus.Rand"] + mdl$VCV[, "units"])
        phylo.df[n, "interMean"] <- mean(inter)
        phylo.df[n, "inter_lowerHDP"] <- HPDinterval(inter)[1] # 95% credible interval 
        phylo.df[n, "inter_upperHDP"] <- HPDinterval(inter)[2]
        ## Intra-genus
        intra <- mdl$VCV[, "units"]/(mdl$VCV[, "genus.Rand"] + mdl$VCV[, "units"])
        phylo.df[n, "intraMean"] <- mean(intra)
        phylo.df[n, "intra_lowerHDP"] <- HPDinterval(intra)[1] # 95% credible interval 
        phylo.df[n, "intra_upperHDP"] <- HPDinterval(intra)[2]
      }
    }
  }
  # Order results
  phylo.df <- arrange(phylo.df, model)
  phylo.df <- unique(phylo.df)
  rslts$phylo_var <-   phylo.df
  return(rslts)
}

# Multi-models indexes
multiModelsIndexFun <- function() {
  rslts <- list()
  phylo.df <- data.frame()
  corr.df <- data.frame()
  phylo.res <- data.frame()
  corr.res <- data.frame()
  # namesvec <- names(multiMod.res)
  # multiMdls <- multiMdls %>% arrange(namesvec)
  for(n in 1:length(multiMdls$type)){ 
    mdl <- multiMod.res[[n]]$"model"
    mdlName <- multiMod.res[[n]]$"name"
    vars <- as.vector(unlist(strsplit(as.character(multiMdls[n, "variables"]), ", ")))[c(1,2)]
  if(is.na(str_extract(mdlName, "NP"))){
    for(v in 1:length(vars)){
      var1 <- vars[v]
      phylo.df[v, "model"] <- mdlName
      phylo.df[v, "variable"] <- var1
      ## Phylogenetic signal
      lambda <- mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")]/
        (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")] + 
           mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
      pd <- p_direction(lambda)
      pd <- as.numeric(pd$pd)
      pval <- 2*(1 - pd)
      phylo.df[v, "lambdaMean"] <- mean(lambda)
      phylo.df[v, "lambda_lowerHPD"] <- HPDinterval(lambda)[1] # 95% credible interval 
      phylo.df[v, "lambda_upperHPD"] <- HPDinterval(lambda)[2]
      phylo.df[v, "lambda_p_value"] <- pval
      ## inter-genus
      inter <- mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")]/
        (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")] + 
           mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
      pd <- p_direction(inter)
      pd <- as.numeric(pd$pd)
      pval <- 2*(1 - pd)
      phylo.df[v, "interMean"] <- mean(inter)
      phylo.df[v, "inter_lowerHDP"] <- HPDinterval(inter)[1] # 95% credible interval 
      phylo.df[v, "inter_upperHDP"] <- HPDinterval(inter)[2]
      phylo.df[v, "inter_p_value"] <- pval
      ## Intra-genus
      intra <- mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")]/
        (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[,paste0("trait", var1, ":trait", var1, ".genus.Rand")] + 
           mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
      pd <- p_direction(intra)
      pd <- as.numeric(pd$pd)
      pval <- 2*(1 - pd)
      phylo.df[v, "intraMean"] <- mean(intra)
      phylo.df[v, "intra_lowerHDP"] <- HPDinterval(intra)[1] # 95% credible interval 
      phylo.df[v, "intra_upperHDP"] <- HPDinterval(intra)[2]
      phylo.df[v, "intra_p_value"] <- pval
      phylo.res <- rbind(phylo.res, phylo.df)
        # Correlated evolution
        for(e in 1:length(vars)){
          var2 <- vars[e]
          corr.df[e, "model"] <- mdlName
          corr.df[e, "var1"] <- var1
          corr.df[e, "var2"] <- var2
          if(var1 == var2){
            next()
          }
          evoCor  <- mdl$VCV[, paste0("trait", var1, ":trait", var2, ".animal")]/
            sqrt(mdl$VCV[, paste0("trait", var2, ":trait", var2, ".animal")] * mdl$VCV[, paste0("trait", var1, ":trait", var1, ".animal")]) 
          # densityplot(evoCor)
          pd <- p_direction(evoCor)
          pd <- as.numeric(pd$pd)
          pval <- 2*(1 - pd)
          corr.df[e, "meanCor"]  <- mean(evoCor)
          corr.df[e, "cor_lowerHDP"]  <- HPDinterval(evoCor, 0.95)[1]
          corr.df[e, "cor_upperHDP"]  <- HPDinterval(evoCor, 0.95)[2]
          corr.df[e, "p_value"] <- pval
          corr.res <- rbind(corr.res, corr.df)
        }
    }
    }
  }
  # Order results
  phylo.res <- arrange(phylo.res, variable, model)
  phylo.res <- unique(phylo.res)
  rslts$phylo_var <-   phylo.res
  corr.res <- corr.res %>% distinct(corr.res$meanCor, .keep_all = T)
  corr.res$`corr.res$meanCor` <- NULL
  corr.res <- arrange(corr.res, var1, var2, model)
  rslts$evoCor <- corr.res
  return(rslts)
}

## Singificant ecoevolutionary correlations
singEcoFun <- function(){
  eCor.res <- read.csv(paste0(output.dir, "EvoCorr_multiModels.csv"), header = T)
  signEco <- data.frame()
  for(i in 1:length(eCor.res$model)){
    mdl <- eCor.res[i, "model"]
    if(isTRUE(eCor.res$cor_lowerHDP[i] < 0 && eCor.res$cor_upperHDP[i] < 0) | isTRUE(eCor.res$cor_lowerHDP[i] > 0 && eCor.res$cor_upperHDP[i] > 0 )){
      signEco <- rbind(signEco, eCor.res[eCor.res$model == mdl, ])
      signEco$model <- str_replace_all(signEco$model, ", ", "_")
    }
  }
  return(signEco)
}


#### UNIRESPONSE MODELS --------------------------------------------------------------------------------------- ####

print("R2 calculation for uni-response models")
# R2s
uniModR2.res <- uniModR2Fun(df.res = uniMod.res)
NP_uniModR2.res <- uniModR2Fun(df.res = NP_uniMod.res)
uniModR2.res
write.csv(uniModR2.res, paste0(output.dir, "R2_uniModels.csv"), row.names = F)
NP_uniModR2.res
write.csv(NP_uniModR2.res, paste0(output.dir, "NP_R2_uniModels.csv"), row.names = F)
print(paste("==>", output.dir, "NP_R2_uniModels.csv"))

# DIC (ordered)
print("DIC calculation for uni-response models")
uniMdls_DIC <- DICFun(uniMdls, uniMod.res)
NP_uniMdls_DIC <- DICFun(uniMdls, NP_uniMod.res)
uniMdls_DIC$NP_DIC <- NP_uniMdls_DIC$DIC
write.csv(uniMdls_DIC, paste0(output.dir, "descriptionAndDIC_uniModels.csv"), row.names = F)
print(paste("==>", output.dir, "description_uniModels.csv"))

# Best model

bestUniMdls <- bestFun(uniMdls_DIC)

#Phylogenetic signal
print("Phylogenetic signal calculation for uni-response models")
uniModelsIndex.res <- uniModelsIndexFun()
uniPhylo.res <- uniModelsIndex.res$phylo_var
write.csv(uniPhylo.res, paste0(output.dir, "phylo_signal_uniModels.csv"), row.names = F)
print(paste("==>", output.dir, "phylo_signal_uniModels.csv"))

# Non-phylogenetic variances
print("Variance components calculation for non-phylogenetic uni-response models")
NP_uniModelsIndex.res <- NP_uniModelsIndexFun()
NP_uniPhylo.res <- NP_uniModelsIndex.res$phylo_var
write.csv(NP_uniPhylo.res, paste0(output.dir, "NP_variances_uniModels.csv"), row.names = F)
print(paste("==>", output.dir, "NP_variances_uniModels.csv"))


#### MULTIRESPONSE MODELS ------------------------------------------------------------------------------------- ####

# DIC
print("DIC calculation for mutli-response models")
multiMdls_DIC <- DICFun(multiMdls, multiMod.res)
# NP_multiMdls_DIC <- DICFun(multiMdls, NP_multiMod.res)
# multiMdls_DIC$NP_DIC <- NP_multiMdls_DIC$DIC
write.csv(multiMdls_DIC, paste0(output.dir, "description_multiModels.csv"), row.names = F)
print(paste("==>", output.dir, "description_multiModels.csv"))

# Best models

bestMultiMdls <- bestFun(mdls = multiMdls_DIC, mdlsType = "Var_covariance_env")
write.csv(bestMultiMdls, paste0(manuscript.dir, "tables/bestMultiModels/bestModelsList.csv"), row.names = F)
print(paste0("===>", manuscript.dir, "tables/bestMultiModels/bestModelsList.csv"))

# Indexes
print("Phylogenetic indexes calculation for multi-response models")
multiModelsIndex.res <- multiModelsIndexFun()
phylo.res <- multiModelsIndex.res$phylo_var
phylo.res
write.csv(phylo.res, paste0(output.dir, "phylo_signal_multiModels.csv"), row.names = F)
print(paste("==>", output.dir, "phylo_signal_multiModels.csv"))

eCor.res <- multiModelsIndex.res$evoCor
eCor.res
write.csv(eCor.res, paste0(output.dir, "EvoCorr_multiModels.csv"), row.names = F)
print(paste("==>", output.dir, "EvoCorr_multiModels.csv"))

signEco <- singEcoFun()
signEco
write.csv(signEco, paste0(output.dir, "SignEvoCorr_multiModels.csv"), row.names = F)
print(paste("==>", output.dir, "SignEvoCorr_multiModels.csv"))
