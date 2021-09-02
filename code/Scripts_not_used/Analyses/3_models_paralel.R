#### MODELS ###################################################################################################################

# P Sanchez-Martinez

source("code/manuscript/1_init.R")

print("Running models in parallel, this script can report errors as it seems that R is not that good doing work parallely...")

#### CORES AND PARALLEL WORK -----------------------------------------------------------------------------------------####

numCores <- detectCores() - 1 
print(paste0("Working with ", numCores, " cores."))
cl <- makeCluster(numCores)
registerDoParallel(cl)

#### FUNCTIONS ---------------------------------------------------------------------------------------------------------####

### Modelling functions

# Uni-response mcmcglmm Function
mcmcglmmFun <- function(dta, fix.frml, ran.frml, fmly, nitte = nitter, burn = burni, thinn = thinni){
        # Model
          # Non-Phylogenetic
          if(is.na(stringr::str_extract(ran.frml, "animal"))){
            prior <- prior1
            mdl <- MCMCglmm::MCMCglmm(as.formula(as.character(fix.frml)), random = as.formula(as.character(ran.frml)),
                            family = fmly, prior = prior1,
                            data = dta, nitt = nitter, burnin = burni, thin = thinni) 
          } else {
          # Phylogenetic
          prior <- prior2
          mdl <- MCMCglmm::MCMCglmm(as.formula(as.character(fix.frml)), random = as.formula(as.character(ran.frml)),
                          family = fmly, pedigree=phylo, prior = prior,
                          data = dta, nitt = nitter, burni = burn, thin = thinni)
          }
  rslts <- list()
  # Model results
  rslts$NumObs <- length(dta[, 1])
  rslts$model <- mdl
  rslts$summary <- summary(mdl)
  # Model diagnostics
  rslts$autocFix <- coda::autocorr.diag(mdl$Sol)
  rslts$autocRan <- coda::autocorr.diag(mdl$VCV)
  rslts$heidelFix <- coda::heidel.diag(mdl$Sol)
  rslts$heidelRan <- coda::heidel.diag(mdl$VCV)
  rslts$effSizeFix <- coda::effectiveSize(mdl$VCV)
  rslts$effSizeRan <- coda::effectiveSize(mdl$Sol)
  # DIC
  rslts$DIC <- mdl$DIC
  return(rslts)
}

# Multi-response mcmcglmm Function
multiMcmcglmmFun <- function(dta, fix.frml, ran.frml, fmly, nitte = nitter, burn = burni, thinn = thinni){
  # Model
  # Non-Phylogenetic
  if(is.na(stringr::str_extract(ran.frml, "animal"))){
    prior <- priorMulti2_NP
    mdl <- MCMCglmm::MCMCglmm(as.formula(as.character(fix.frml)), 
                    random = as.formula(as.character(ran.frml)), 
                    rcov = ~us(trait):units, family = fmly,
                    prior = prior, data= dta, nitt = nitter, burnin = burni, thin = thinni)
  } else {
    # Phylogenetic
    prior <- priorMulti2
    mdl <- MCMCglmm::MCMCglmm(as.formula(as.character(fix.frml)), 
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
  rslts$autocFix <- coda::autocorr.diag(mdl$Sol)
  rslts$autocRan <- coda::autocorr.diag(mdl$VCV)
  rslts$heidelFix <- coda::heidel.diag(mdl$Sol)
  rslts$heidelRan <- coda::heidel.diag(mdl$VCV)
  rslts$effSizeFix <- coda::effectiveSize(mdl$VCV)
  rslts$effSizeRan <- coda::effectiveSize(mdl$Sol)
  # DIC
  rslts$DIC <- mdl$DIC
  return(rslts)
}


# Uni-response models
uniResponseModels <- function(mdls, randForm){
  uniMod.res <- list()
  uniMod.res <- foreach(i = 1:length(mdls$type), .export=ls(envir=globalenv())) %dopar% {
    mdl <- paste0(as.character(mdls[i, "type"]), "_", as.character(mdls[i, "group"]), "_", as.character(mdls[i, "variables"]))
    mdl <- stringr::str_replace(mdl, ", ", "_")
    if(randForm == "NP_ran.frml"){
      mdl <- paste0("NP_", mdl)
    }
    # Complete data
    vars <- as.vector(unlist(strsplit(as.character(mdls[i, "variables"]), ", ")))
    completeData <- completePhyloFun(plmm.data, desiredCols = vars)
    # Model
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    if(mdls[i, "N_ResVars"] == 1){ 
      fmly <- fmly 
    }
    mcmcglmmFun(dta = completeData, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
  }
  return(uniMod.res)
}

  
# Multi-response models
multiResponseModels <- function(mdls, randForm){
  multiMod.res <- list()
  multiMod.res <- foreach(i = 1:length(mdls$type), .export=ls(envir=globalenv())) %dopar% {
    mdl <- paste0(as.character(mdls[i, "type"]), "_", as.character(mdls[i, "group"]), "_", as.character(mdls[i, "variables"]))
    mdl <- stringr::str_replace(mdl, ", ", "_")
    if(randForm == "NP_ran.frml"){
    mdl <- paste0("NP_", mdl)
      }
    # Complete data
    vars <- as.vector(unlist(strsplit(as.character(mdls[i, "variables"]), ", ")))
    completeData <- completePhyloFun(plmm.data, desiredCols = vars)
    # Model 
    fix.frml <- mdls[i, "fix.frml"]
    ran.frml <- mdls[i, randForm]
    fmly <- fmly2
    multiMcmcglmmFun(dta = completeData, fix.frml = fix.frml, ran.frml = ran.frml, fmly = fmly, nitte = nitter, burn = burni, thinn = thinni)
  }
  # Model names
  for(i in 1:length(mdls$type)){
  mdl <- paste0(as.character(mdls[i, "type"]), "_", as.character(mdls[i, "group"]), "_", as.character(mdls[i, "variables"]))
  mdl <- stringr::str_replace(mdl, ", ", "_")
  if(randForm == "NP_ran.frml"){
    mdl <- paste0("NP_", mdl)
  }
  multiMod.res[[i]]$name <- mdl
  }
  return(multiMod.res)
}


#### UNI-RESPONSE MODELS ----------------------------------------------------------------------------------------------####

tic()

### Phylogenetic 
uniMod.res <- uniResponseModels(uniMdls, randForm = "ran.frml")

### Non-phylogenetic
# NP_uniMod.res <- uniResponseModels(uniMdls, randForm = "NP_ran.frml")


#### MULTI-RESPONSE MODELS --------------------------------------------------------------------------------------------####

## Phylogenetic
multiMod.res <- multiResponseModels(multiMdls, randForm = "ran.frml")

# ## Non phylogenetic
NP_multiMod.res <- multiResponseModels(multiMdls, randForm = "NP_ran.frml")

# # Save image
# save.image(file = "data/RData/plmm_results.RData")
# print("====> data/RData/plmm_results.RData")
toc()