#### MODELS ###################################################################################################################

# Pablo Sanchez
# 10-2019

setwd("/home/jordi/Pablo")
source("code/manuscript/1_init.R")

print("Running models, the whole script may take few weeks of computation, depending on models specifications...")

#### FUNCTIONS ---------------------------------------------------------------------------------------------------------####

### Auxiliary functions

# Complete data
completePhyloFun <- function(dataset, desiredCols) {
  completeVec <- complete.cases(dataset[, desiredCols])
  completeData <- dataset[completeVec, ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% completeData$genus])
  completeData <- completeData[completeData$genus %in% phylo$tip.label, ]
  return(as.data.frame(completeData))
}

### Modelling functions

# Uni-response mcmcglmm Function
mcmcglmmFun <- function(dta, fix.frml, ran.frml, fmly, nitte = nitter, burn = burni, thinn = thinni){
        # Model
          # Non-Phylogenetic
          if(is.na(str_extract(ran.frml, "animal"))){
            prior <- prior1
            mdl <- MCMCglmm(as.formula(as.character(fix.frml)), random = as.formula(as.character(ran.frml)),
                            family = fmly, prior = prior1,
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

# Multi-response mcmcglmm Function
multiMcmcglmmFun <- function(dta, fix.frml, ran.frml, fmly, nitte = nitter, burn = burni, thinn = thinni){
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

# Uni-response models
uniResponseModels <- function(mdls, randForm){
  uniMod.res <- list()
  for(i in 1:length(mdls$type)){
    mdl <- paste0(as.character(mdls[i, "type"]), "_", as.character(mdls[i, "group"]), "_", as.character(mdls[i, "variables"]))
    mdl <- str_replace(mdl, ", ", "_")
    if(randForm == "NP_ran.frml"){
      mdl <- paste0("NP_", mdl)
    }
    print(mdl)
    # Complete data
    vars <- as.vector(unlist(strsplit(as.character(mdls[i, "variables"]), ", ")))
    completeData <- completePhyloFun(plmm.data, desiredCols = vars)
    # Model 
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    if(mdls[i, "N_ResVars"] == 1){ 
      fmly <- fmly 
      }
    uniMod.res[[mdl]] <- mcmcglmmFun(dta = completeData, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
    uniMod.res[[mdl]]$"name" <- mdl
  }
  return(uniMod.res)
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
    completeData <- completePhyloFun(plmm.data, desiredCols = vars)
    # Model 
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    fmly <- fmly2
    multiMod.res[[mdl]] <- multiMcmcglmmFun(dta = completeData, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
    multiMod.res[[mdl]]$"name" <- mdl
  }
  return(multiMod.res)
}


#### UNI-RESPONSE MODELS ----------------------------------------------------------------------------------------------####

tic()

### Phylogenetic 
# uniMod.res <- uniResponseModels(uniMdls, randForm = "ran.frml")

### Non-phylogenetic
# NP_uniMod.res <- uniResponseModels(uniMdls, randForm = "NP_ran.frml")


#### MULTI-RESPONSE MODELS --------------------------------------------------------------------------------------------####

## Phylogenetic
multiMod.res <- multiResponseModels(multiMdls, randForm = "ran.frml")

## Non phylogenetic
NP_multiMod.res <- multiResponseModels(multiMdls, randForm = "NP_ran.frml")

# Save image
save.image(file = "data/RData/plmm_results.RData")
print("====> data/RData/plmm_results.RData")
toc()