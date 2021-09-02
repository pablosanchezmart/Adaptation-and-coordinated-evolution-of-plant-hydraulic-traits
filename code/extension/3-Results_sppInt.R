#### MODELS RESULTS ###################################################################################################################

# P Sanchez-Martinez

print("Running results checking for potential effects of species number...")
remove(list = ls())

print("loading multivariate models results...")
load("data/RData/plmm_multivariate_results_sppInt.RData")

source("code/extension/sppInt_init.R")

### FUNCTIONS ###

# DIC
DICFun <- function(mdls, mod.res){
  mdls$DIC <- NA
  for(i in 1:length(mdls$type)){
    mdls$DIC[i] <- mod.res[[i]]$DIC
  }
  mdls <- arrange(mdls,variables, desc(DIC))
  return(mdls)
}

# Phylogenetic indexes

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
singEcoFun <- function(eCor.res){
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


#### MULTIRESPONSE MODELS ------------------------------------------------------------------------------------- ####

# DIC
print("DIC calculation for mutli-response models")
multiMdls_DIC <- DICFun(multiMdls, multiMod.res)
# NP_multiMdls_DIC <- DICFun(multiMdls, NP_multiMod.res)
# multiMdls_DIC$NP_DIC <- NP_multiMdls_DIC$DIC

# Indexes
print("Phylogenetic indexes calculation for multi-response models")
multiModelsIndex.res <- multiModelsIndexFun()
phylo.res <- multiModelsIndex.res$phylo_var
phylo.res

eCor.res <- multiModelsIndex.res$evoCor
eCor.res
write.csv(eCor.res, paste0(output.dir, "EvoCorr_multiModels_sppInt.csv"), row.names = F)
print(paste("==>", output.dir, "EvoCorr_multiModels_sppInt.csv"))

signEco <- singEcoFun(eCor.res)
signEco
write.csv(signEco, paste0(output.dir, "SignEvoCorr_multiModels_sppInt.csv"), row.names = F)
print(paste("==>", output.dir, "SignEvoCorr_multiModels_sppInt.csv"))
