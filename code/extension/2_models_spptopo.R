#### MODELS ###################################################################################################################

# P Sanchez-Martinez

print("Running models, the whole script may take few weeks of computation, depending on models specifications...")

#### FUNCTIONS ---------------------------------------------------------------------------------------------------------####

### Modelling functions

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
    completeData <- completePhyloFun(tree = ht.tree, dataset = plmm.data, desiredCols = vars)
    # Model 
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    fmly <- fmly2
    multiMod.res[[mdl]] <- multiMcmcglmmFun(dta = completeData$dta, phylo = completeData$phylo, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
    multiMod.res[[mdl]]$"name" <- mdl
  }
  return(multiMod.res)
}


#### MULTI-RESPONSE MODELS --------------------------------------------------------------------------------------------####

## Phylogenetic
source("code/extension/spptopo_init.R")
tic()
multiMod.res <- multiResponseModels(multiMdls, randForm = "ran.frml")

# Save image
# save.image(file = "data/RData/plmm_multivariate_results_spptopo.RData")
# print("====> data/RData/plmm_multivariate_results_spptopo.RData")
# toc()

remove(list = ls())
