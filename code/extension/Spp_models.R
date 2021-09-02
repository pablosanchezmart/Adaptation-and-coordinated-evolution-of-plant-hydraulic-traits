#### MODELS ###################################################################################################################

# P Sanchez-Martinez

print("Running models, the whole script may take few weeks of computation, depending on models specifications...")

source("code/extension/Spp_init.R")

#### FUNCTIONS ---------------------------------------------------------------------------------------------------------####

### Modelling functions

# Uni-response mcmcglmm Function
uniMcmcglmmFun <- function(dta, phylo, fix.frml, ran.frml, fmly, nitte = nitter, burn = burni, thinn = thinni){
  # Model
  # Non-Phylogenetic
  if(is.na(str_extract(ran.frml, "animal"))){
    prior <- prior2_NP
    mdl <- MCMCglmm(as.formula(as.character(fix.frml)),
                    family = fmly, prior = prior,
                    data = dta, nitt = nitter, burnin = burni, thin = thinni) 
  } else {
    # Phylogenetic
    prior <- prior2
    mdl <- MCMCglmm(as.formula(as.character(fix.frml)), random = as.formula(as.character(ran.frml)),
                    family = fmly, pedigree=phylo, prior = prior,
                    data = dta, nitt = nitter, burni = burn, thin = thinni)
  }
  rslts <- list()
  # Model results
  rslts$NumObs <- length(dta[, 1])
  rslts$model <- mdl
  rslts$summary <- summary(mdl)
  # Model diagnostics
  rslts$autocFix <- autocorr.diag(mdl$Sol)
  rslts$autocRan <- autocorr.diag(mdl$VCV)
  rslts$heidelFix <- heidel.diag(mdl$Sol)
  rslts$heidelRan <- heidel.diag(mdl$VCV)
  rslts$effSizeFix <- effectiveSize(mdl$VCV)
  rslts$effSizeRan <- effectiveSize(mdl$Sol)
  # DIC
  rslts$DIC <- mdl$DIC
  return(rslts)
}

# Uni-response models
uniResponseModels <- function(mdls, randForm){
  multiMod.res <- list()
  for(i in 1:length(mdls$type)){
    mdl <- paste0(as.character(mdls[i, "type"]), "_", as.character(mdls[i, "group"]), "_", as.character(mdls[i, "variables"]))
    mdl <- str_replace(mdl, ", ", "_")
    if(randForm == "NP_ran.frml"){
      mdl <- paste0("NP_", mdl)
    }
    print(mdl)
    # Complete data
    vars <- as.vector(unlist(strsplit(as.character(mdls[i, "variables"]), ", ")))
    if(!is.na(str_extract(mdls[i, "fix.frml"], "group_leaf"))){
      completeData <- completePhyloFun(tree = ht.tree, dataset = plmm.data, desiredCols = c(vars, "group_leaf"))
    } else {
      completeData <- completePhyloFun(tree = ht.tree, dataset = plmm.data, desiredCols = vars) 
    }
    # Model 
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    fmly <- fmly
    multiMod.res[[mdl]] <- uniMcmcglmmFun(dta = completeData$dta, phylo = completeData$phylo, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
    multiMod.res[[mdl]]$"name" <- mdl
  }
  return(multiMod.res)
}

# Multi-response mcmcglmm Function
multiMcmcglmmFun <- function(dta, phylo, fix.frml, ran.frml, fmly, nitte = nitter, burn = burni, thinn = thinni){
  # Model
  # Non-Phylogenetic
  if(is.na(str_extract(ran.frml, "animal"))){
    prior <- priorMulti2_NP
    mdl <- MCMCglmm(as.formula(as.character(fix.frml)), 
                    random = as.formula(as.character(ran.frml)), 
                    rcov = ~us(trait):units, family = fmly,
                    prior = prior, data= dta, nitt = nitter, burnin = burni, thin = thinni)
  } else {
    # Phylogenetic
    prior <- priorMulti2
    mdl <- MCMCglmm(as.formula(as.character(fix.frml)), 
                    random = as.formula(as.character(ran.frml)), 
                    rcov = ~us(trait):units, pedigree=phylo, family = fmly,
                    prior = prior, data= dta, nitt = nitter, burnin = burni, thin = thinni)
  }
  rslts <- list()
  # Model results
  rslts$NumObs <- length(dta[, 1])
  rslts$model <- mdl
  rslts$summary <- summary(mdl)
  # Model diagnostics
  rslts$autocFix <- autocorr.diag(mdl$Sol)
  rslts$autocRan <- autocorr.diag(mdl$VCV)
  rslts$heidelFix <- heidel.diag(mdl$Sol)
  rslts$heidelRan <- heidel.diag(mdl$VCV)
  rslts$effSizeFix <- effectiveSize(mdl$VCV)
  rslts$effSizeRan <- effectiveSize(mdl$Sol)
  # DIC
  rslts$DIC <- mdl$DIC
  return(rslts)
}

# Multi-response models
multiResponseModels <- function(mdls, randForm){
  multiMod.res <- list()
  for(i in 1:length(mdls$type)){
    mdl <- paste0(as.character(mdls[i, "type"]), "_", as.character(mdls[i, "group"]), "_", as.character(mdls[i, "variables"]))
    mdl <- str_replace(mdl, ", ", "_")
    if(randForm == "NP_ran.frml"){
      mdl <- paste0("NP_", mdl)
    }
    print(mdl)
    # Complete data
    vars <- as.vector(unlist(strsplit(as.character(mdls[i, "variables"]), ", ")))
    if(!is.na(str_extract(mdls[i, "fix.frml"], "group_leaf"))){
      completeData <- completePhyloFun(tree = ht.tree, dataset = plmm.data, desiredCols = c(vars, "group_leaf"))
    } else {
      completeData <- completePhyloFun(tree = ht.tree, dataset = plmm.data, desiredCols = vars) 
    }
    # Model 
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    fmly <- fmly2
    multiMod.res[[mdl]] <- multiMcmcglmmFun(dta = completeData$dta, phylo = completeData$phylo, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
    multiMod.res[[mdl]]$"name" <- mdl
  }
  return(multiMod.res)
}


#### UNI-RESPONSE MODELS --------------------------------------------------------------------------------------------####

## Phylogenetic
# Spp_uniMod.res <- uniResponseModels(mdls = Sp_uniMdls, randForm = "ran.frml")
# 
# # Save image
# save.image(file = "data/RData/plmm_univariate_Spp_results.RData")
# print("====> data/RData/plmm_univariate_Spp_results.RData.RData")
# toc()

## Non-Phylogenetic
# Spp_uniMod_NP.res <- uniResponseModels(mdls = Sp_uniMdls, randForm = "NP_ran.frml")
# 
# # Save image
# save.image(file = "data/RData/plmm_univariate_NP_Spp_results.RData")
# print("====> data/RData/plmm_univariate_NP_Spp_results.RData")
# toc()


#### MULTI-RESPONSE MODELS --------------------------------------------------------------------------------------------####

# ## Phylogenetic
# Spp_multiMod.res <- multiResponseModels(mdls = Sp_multiMdls, randForm = "ran.frml")
# 
# # Save image
# save.image(file = "data/RData/plmm_multivariate_Spp_results.RData")
# print("====> data/RData/plmm_multivariate_Spp_results.RData")
# toc()

