#### SPECIES LEVEL PHYLOGENY MODELS RESULTS ###################################################################################################################

# P Sanchez-Martinez

print("Running Spp level results...")
remove(list = ls())

print("loading univariate models results...")
load("data/RData/plmm_univariate_Spp_results_diagn.RData")

print("loading multivariate models results...")
load("data/RData/plmm_multivariate_Spp_results_diagn.RData")

source("code/extension/Spp_init.R")

### FUNCTIONS ###

# Phylogenetic indexes

# Uni-models phylogenetic signal
Spp_uniModelsIndexFun <- function() {
  rslts <- list()
  phylo.df <- data.frame()
  for(n in 1:length(Sp_uniMdls$type)){ 
    var1 <- as.vector(unlist(strsplit(as.character(Sp_uniMdls[n, "variables"]), ", ")))[1]
    mdlName <- Sp_uniMdls[n, "name"]
    mdl <- Spp_uniMod.res[[mdlName]]$model
    nObs <- Spp_uniMod.res[[mdlName]]$NumObs
    if(is.na(str_extract(mdlName, "NP"))){
      phylo.df[n, "model"] <- mdlName
      phylo.df[n, "variable"] <- var1
      phylo.df[n, "N"] <- nObs
      ## Phylogenetic signal
      lambda <- mdl$VCV[,"animal"]/ (mdl$VCV[, "animal"] + mdl$VCV[, "units"])
      phylo.df[n, "lambdaMean"] <- mean(lambda)
      phylo.df[n, "lambda_lowerHPD"] <- HPDinterval(lambda)[1] # 95% credible interval 
      phylo.df[n, "lambda_upperHPD"] <- HPDinterval(lambda)[2]
      ## Intra-species
      intra <- mdl$VCV[, "units"]/(mdl$VCV[, "animal"] + mdl$VCV[, "units"])
      phylo.df[n, "intraMean"] <- mean(intra)
      phylo.df[n, "intra_lowerHDP"] <- HPDinterval(intra)[1] # 95% credible interval 
      phylo.df[n, "intra_upperHDP"] <- HPDinterval(intra)[2]
    }
  }
  # Order results
  phylo.df <- arrange(phylo.df, model)
  phylo.df <- unique(phylo.df)
  rslts$phylo_var <-   phylo.df
  return(rslts)
}

# Multi-models indexes
Spp_multiModelsIndexFun <- function() {
  rslts <- list()
  phylo.df <- data.frame()
  corr.df <- data.frame()
  phylo.res <- data.frame()
  corr.res <- data.frame()
  namesvec <- names(Spp_multiMod.res)
  for(n in 1:length(namesvec)){ 
    mdl <- Spp_multiMod.res[[n]]$"model"
    mdlName <- Spp_multiMod.res[[n]]$"name"
    vars <- as.vector(unlist(strsplit(as.character(Sp_multiMdls[Sp_multiMdls$name == mdlName, "variables"]), ", ")))[c(1,2)]
    if(is.na(str_extract(mdlName, "NP"))){
      for(v in 1:length(vars)){
        var1 <- vars[v]
        phylo.df[v, "model"] <- mdlName
        phylo.df[v, "variable"] <- var1
        ## Phylogenetic signal
        lambda <- mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")]/
          (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
        pd <- p_direction(lambda)
        pd <- as.numeric(pd$pd)
        pval <- 2*(1 - pd)
        phylo.df[v, "lambdaMean"] <- mean(lambda)
        phylo.df[v, "lambda_lowerHPD"] <- HPDinterval(lambda)[1] # 95% credible interval 
        phylo.df[v, "lambda_upperHPD"] <- HPDinterval(lambda)[2]
        phylo.df[v, "p_direction"] <- pd * 100
        phylo.df[v, "p_value"] <- pval
        ## Intra-species
        intra <- mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")]/
          (mdl$VCV[,paste0("trait", var1, ":trait", var1, ".animal")] + mdl$VCV[, paste0("trait", var1, ":trait", var1, ".units")])
        phylo.df[v, "intraMean"] <- mean(intra)
        phylo.df[v, "intra_lowerHDP"] <- HPDinterval(intra)[1] # 95% credible interval 
        phylo.df[v, "intra_upperHDP"] <- HPDinterval(intra)[2]
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
          pd <- p_direction(evoCor)
          pd <- as.numeric(pd$pd)
          pval <- 2*(1 - pd)
          corr.df[e, "meanCor"]  <- mean(evoCor)
          corr.df[e, "cor_lowerHDP"]  <- HPDinterval(evoCor, 0.95)[1]
          corr.df[e, "cor_upperHDP"]  <- HPDinterval(evoCor, 0.95)[2]
          corr.df[e, "p_direction"] <- pd * 100
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

#### UNIRESPONSE MODELS --------------------------------------------------------------------------------------- ####

#Phylogenetic signal
print("Phylogenetic signal calculation for uni-response models")
uniModelsIndex.res <- Spp_uniModelsIndexFun()
uniPhylo.res <- uniModelsIndex.res$phylo_var
write.csv(uniPhylo.res, paste0(output.dir, "extension/Spp_phylo_signal_uniModels.csv"), row.names = F)
print(paste("==>", output.dir, "extension/Spp_phylo_signal_uniModels.csv"))

#### MULTIRESPONSE MODELS ------------------------------------------------------------------------------------- ####

# Indexes
print("Phylogenetic indexes calculation for multi-response models")
multiModelsIndex.res <- Spp_multiModelsIndexFun()
phylo.res <- multiModelsIndex.res$phylo_var
phylo.res
write.csv(phylo.res, paste0(manuscript.dir, "tables/SA2/Spp_phylo_signal_multiModels.csv"), row.names = F)
print(paste("==>", output.dir, "manuscript/tables/SA2/Spp_phylo_signal_multiModels.csv"))

eCor.res <- multiModelsIndex.res$evoCor
eCor.res$model <- str_replace_all(eCor.res$model, ", ", "_")
eCor.res
write.csv(eCor.res, paste0(manuscript.dir, "tables/SA2/Spp_EvoCorr_multiModels.csv"), row.names = F)
print(paste("==>", output.dir, "manuscript/tables/SA2/Spp_EvoCorr_multiModels.csv"))

