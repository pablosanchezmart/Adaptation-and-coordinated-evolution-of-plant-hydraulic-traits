#### MODELS DIAGNOSTICS AND RERUNNING ###################################################################################################################

# P Sanchez-Martinez

print("Running model diagnostics...")

#### UNI-RESPONSE MODELS ----------------------------------------------------------------------------------------------####

### Phylogenetic

source("code/manuscript/1_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_univariate_results_diagn.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
uniMod.diag <- diagnFun(uniMod.res)
uniMod.diag # 2
uniMod.fail <- failedMdls(uniMod.diag)

# Rerun
uniMod.fail <- uniMdls[uniMdls$name %in% uniMod.fail$model, ]
uniMod.fail$model <- uniMod.fail$name

for(mdl in uniMod.fail$model){
  cat("Re-Running ", mdl)
  uniMdl <- uniMdls[which(uniMdls$name == mdl), ]
  uniMod.res[mdl] <- uniResponseModels(uniMdl, randForm = "ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_univariate_results_diagn.RData")
print("====> data/RData/plmm_univariate_results_diagn.RData")
toc()

remove(list = ls())


### Non-phylogenetic

source("code/manuscript/1_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_univariate_NP_results_diagn.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
NP_uniMod.diag <- diagnFun(NP_uniMod.res)
NP_uniMod.diag #1
NP_uniMod.fail <- failedMdls(NP_uniMod.diag)

# Rerun
NP_uniMod.fail <- uniMdls[uniMdls$name %in% str_remove_all(NP_uniMod.fail$model, "NP_"), ]
NP_uniMod.fail$model <- NP_uniMod.fail$name

for(mdl in NP_uniMod.fail$model){
  cat("Re-Running NP_", mdl)
  uniMdl <- uniMdls[which(uniMdls$name == mdl), ]
  mdl.name <- paste0("NP_", mdl)
  NP_uniMod.res[mdl.name] <- uniResponseModels(uniMdl, randForm = "NP_ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_univariate_NP_results_diagn.RData")
print("====> data/RData/plmm_univariate_NP_results_diagn.RData")
toc()

remove(list = ls())

#### MULTI-RESPONSE MODELS --------------------------------------------------------------------------------------------####

## Phylogenetic

source("code/manuscript/1_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_results_diagn.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
multiMod.diag <- diagnFun(multiMod.res)
multiMod.diag #4
multiMod.fail <- failedMdls(multiMod.diag)

# Rerun
multiMod.fail <- multiMdls[multiMdls$name %in% multiMod.fail$model, ]
multiMod.fail$model <- multiMod.fail$name

multiMod.fail <- multiMdls[40:41, ]
multiMod.fail$model <- multiMod.fail$name

for(mdl in multiMod.fail$model){
  cat("Re-Running ", mdl)
  multiMdl <- multiMdls[which(multiMdls$name == mdl), ]
  multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_results_diagn.RData")
print("====> data/RData/plmm_multivariate_results_diagn.RData")
toc()

remove(list = ls())


### Non phylogenetic

source("code/manuscript/1_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_NP_results.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
NP_multiMod.diag <- diagnFun(NP_multiMod.res)
NP_multiMod.diag #7
NP_multiMod.fail <- failedMdls(NP_multiMod.diag)

for(mdl in NP_multiMod.fail$model){
  if(!is.na(str_extract(mdl, "NP"))){
    mdl1 <- str_remove(mdl, "NP_")
  }
  cat("Re-Running ", mdl)
  multiMdl <- multiMdls[which(multiMdls$name == mdl1), ]
  NP_multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "NP_ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_NP_results_diagn.RData")
print("====> data/RData/plmm_multivariate_NP_results_diagn.RData")
toc()

remove(list = ls())


#### UNI-RESPONSE SPECIES LEVEL MODELS (extension) --------------------------------------------------------------------------------------------####

### Phylogenetic

source("code/extension/Spp_init.R")

# load("data/plmm_univariate_Spp_results1.RData")
# Spp_uniMod.res1 <- Spp_uniMod.res
# rm(Spp_uniMod.res)
# load("data/plmm_univariate_Spp_results2.RData")
# Spp_uniMod.res2 <- Spp_uniMod.res
# rm(Spp_uniMod.res)
# Spp_uniMod.res <- c(Spp_uniMod.res1, Spp_uniMod.res2)
# rm(Spp_uniMod.res1)
# rm(Spp_uniMod.res2)
# 
# save.image(file = "data/RData/plmm_univariate_Spp_results.RData")
# print("====> data/RData/plmm_univariate_Spp_results.RData")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("Data/RData/plmm_univariate_Spp_results.RData")

nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
uniMod.diag <- diagnFun(Spp_uniMod.res)
uniMod.diag # 2
uniMod.fail <- failedMdls(uniMod.diag)

# Rerun
uniMod.fail <- Sp_uniMdls[Sp_uniMdls$name %in% uniMod.fail$model, ]
uniMod.fail$model <- uniMod.fail$name

for(mdl in uniMod.fail$model){
  cat("Re-Running ", mdl)
  uniMdl <- Sp_uniMdls[which(Sp_uniMdls$name == mdl), ]
  Spp_uniMod.res[mdl] <- uniResponseModels(uniMdl, randForm = "ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_univariate_Spp_results_diagn.RData")
print("====> data/RData/plmm_univariate_results_diagn.RData")
toc()

remove(list = ls())

## Non-Phylogenetic

source("code/extension/Spp_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_univariate_NP_Spp_results_diagn.RData")

nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
NP_uniMod.diag <- diagnFun(Spp_uniMod_NP.res)
NP_uniMod.diag #4
NP_uniMod.fail <- failedMdls(NP_uniMod.diag)

# Rerun
NP_uniMod.fail <- Sp_uniMdls[Sp_uniMdls$name %in% str_remove_all(NP_uniMod.fail$model, "NP_"), ]
NP_uniMod.fail$model <- NP_uniMod.fail$name

for(mdl in NP_uniMod.fail$model){
  cat("Re-Running NP_", mdl)
  uniMdl <- Sp_uniMdls[which(Sp_uniMdls$name == mdl), ]
  mdl.name <- paste0("NP_", mdl)
  Spp_uniMod_NP.res[mdl.name] <- uniResponseModels(uniMdl, randForm = "NP_ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_univariate_NP_Spp_results_diagn.RData")
print("====> data/RData/plmm_multivariate_NP_Spp_results_diagn.RData")
toc()

remove(list = ls())


#### MULTI-RESPONSE SPECIES LEVEL MODELS (extension) --------------------------------------------------------------------------------------------####

load("data/plmm_multivariate_Spp_results_diagn1.RData")

load("data/plmm_multivariate_Spp_results_diagn2.RData")

load("data/RData/plmm_multivariate_Spp_results_diagn.RData")

Spp_multiMod.res <- c(Spp_multiMod.res, Spp_multiMod.res1, Spp_multiMod.res2)
rm(Spp_multiMod.res1)
rm(Spp_multiMod.res2)

save.image(file = "data/RData/plmm_multivariate_Spp_results_diagn.RData")
print("====> data/RData/plmm_multivariate_Spp_results_diagn.RData")

load("Data/RData/plmm_univariate_Spp_results_diagn.RData")

## Phylogenetic

source("code/extension/Spp_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_Spp_results.RData")

nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
multiMod.diag <- diagnFun(Spp_multiMod.res)
multiMod.diag #4
multiMod.fail <- failedMdls(multiMod.diag)

# Rerun
multiMod.fail <- Sp_multiMdls[Sp_multiMdls$name %in% multiMod.fail$model, ]
multiMod.fail$model <- multiMod.fail$name

for(mdl in multiMod.fail$model){
  cat("Re-Running ", mdl)
  multiMdl <- Sp_multiMdls[which(Sp_multiMdls$name == mdl), ]
  Spp_multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_Spp_results_diagn.RData")
print("====> data/RData/plmm_multivariate_Spp_results_diagn.RData")
toc()

remove(list = ls())




#### MULTI-RESPONSE MODELS (same species that spp lvl phylogeny) --------------------------------------------------------------------------------------------####

## Phylogenetic

source("code/manuscript/1_init.R")

# load("data/plmm_multivariate_results_sppInt1.RData")
# 
# load("data/plmm_multivariate_results_sppInt2.RData")
# 
# multiMod.res <- c(multiMod.res1, multiMod.res2)
# rm(multiMod.res1)
# rm(multiMod.res2)
# 
# save.image(file = "data/RData/plmm_multivariate_results_sppInt.RData")
# print("====> data/RData/plmm_multivariate_results_sppInt.RData")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_results_sppInt.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
multiMod.diag <- diagnFun(multiMod.res)
multiMod.diag #4
multiMod.fail <- failedMdls(multiMod.diag)

# Rerun
multiMod.fail <- multiMdls[multiMdls$name %in% multiMod.fail$model, ]
multiMod.fail$model <- multiMod.fail$name

multiMod.fail <- multiMdls[40:41, ]
multiMod.fail$model <- multiMod.fail$name

for(mdl in multiMod.fail$model){
  cat("Re-Running ", mdl)
  multiMdl <- multiMdls[which(multiMdls$name == mdl), ]
  multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_results_diagn.RData")
print("====> data/RData/plmm_multivariate_results_diagn.RData")
toc()

remove(list = ls())


### Non phylogenetic

source("code/manuscript/1_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_NP_results.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
NP_multiMod.diag <- diagnFun(NP_multiMod.res)
NP_multiMod.diag #7
NP_multiMod.fail <- failedMdls(NP_multiMod.diag)

for(mdl in NP_multiMod.fail$model){
  if(!is.na(str_extract(mdl, "NP"))){
    mdl1 <- str_remove(mdl, "NP_")
  }
  cat("Re-Running ", mdl)
  multiMdl <- multiMdls[which(multiMdls$name == mdl1), ]
  NP_multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "NP_ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_NP_results_diagn.RData")
print("====> data/RData/plmm_multivariate_NP_results_diagn.RData")
toc()

remove(list = ls())


#### MULTI-RESPONSE MODELS (spp level but at genus level, to check for topology effects) ####

## Phylogenetic

source("code/extension/spptopo_init.R")

load("data/plmm_multivariate_results_spptopo1.RData")

load("data/plmm_multivariate_results_spptopo2.RData")

multiMod.res <- c(multiMod.res1, multiMod.res2)
rm(multiMod.res1)
rm(multiMod.res2)

save.image(file = "data/RData/plmm_multivariate_results_sppTopo.RData")
print("====> data/RData/plmm_multivariate_results_sppTopo.RData")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_results_sppTopo.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
multiMod.diag <- diagnFun(multiMod.res)
multiMod.diag #4
multiMod.fail <- failedMdls(multiMod.diag)

# Rerun
multiMod.fail <- multiMdls[multiMdls$name %in% multiMod.fail$model, ]
multiMod.fail$model <- multiMod.fail$name

multiMod.fail <- multiMdls[40:41, ]
multiMod.fail$model <- multiMod.fail$name

for(mdl in multiMod.fail$model){
  cat("Re-Running ", mdl)
  multiMdl <- multiMdls[which(multiMdls$name == mdl), ]
  multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_results_diagn.RData")
print("====> data/RData/plmm_multivariate_results_diagn.RData")
toc()

remove(list = ls())


### Non phylogenetic

source("code/manuscript/1_init.R")

### Functions

## Model diagnostics
diagnFun <- function(mdls.res){
  mdlDiagn <- data.frame("model" = rep(NA, length(mdls.res)),"AutcorrFix" = NA, "AutcorrRan" = NA, "HeidFix" = NA, "heidRan" = NA, "effSizeFix" = NA, "effSizeRan" = NA)
  for(n in 1:length(mdls.res)){ 
    mdl <- mdls.res[[n]]
    name <- mdls.res[[n]]$"name"
    mdlDiagn$model[n] <- name
    # Autocorrelation
    if(any(mdl$autocFix[-1, ] > 0.1)){
      cat(name, " autocFix Failed \n")
      mdlDiagn$AutcorrFix[n] <- "F"
    } else {
      mdlDiagn$AutcorrFix[n] <- "Passed"
    }
    if(any(mdl$autocRan[-1, ] > 0.1)){
      cat(name, " autocRan Failed \n")
      mdlDiagn$AutcorrRan[n] <- "F"
    } else {
      mdlDiagn$AutcorrRan[n] <- "Passed"
    }
    # Convergence (stationary)
    i <- length(mdl$summary$solutions)/5
    heidelFix <- mdl$heidelFix[1:i, 1:3]
    if(any(heidelFix == 0)){
      cat(name, " heidelFix Failed \n")
      mdlDiagn$HeidFix[n] <- "F"
    } else {
      mdlDiagn$HeidFix[n] <- "Passed"
    }
    i <- (length(mdl$summary$Gcovariances)/4) + 1
    heidelRan <- mdl$heidelRan[1:i, 1:3]
    if(any(heidelRan == 0)){
      cat(name, " heidelRan Failed \n")
      mdlDiagn$heidRan[n] <- "F"
    } else {
      mdlDiagn$heidRan[n] <- "Passed"
    }
    # Effect size
    if(any(mdl$effSizeFix < 1000)){
      cat(name, " EffSizeFix Failed \n")
      mdlDiagn$effSizeFix[n] <- "F"
    } else {
      mdlDiagn$effSizeFix[n] <- "Passed"
    }
    if(any(mdl$effSizeRan < 1000)){
      cat(name, " EffSizeRan Failed \n")
      mdlDiagn$effSizeRan[n] <- "F"
    } else {
      mdlDiagn$effSizeRan[n] <- "Passed"
    }
  }
  return(mdlDiagn)
}

## Failed models to re-run increasing sample size
failedMdls <- function(diag.df){
  failed.df <- data.frame()
  for(i in 1:length(diag.df[, "model"])){
    if(any(diag.df[i, ] == "F")){
      failed.df <- rbind(failed.df, diag.df[i, ])
    }
  }
  return(failed.df)
}

# Rdata
tic()
load("data/RData/plmm_multivariate_NP_results.RData")

nitter = 15000000
burni = 1000000
thinni = 4000
set.seed(123)

# Diagnostics
NP_multiMod.diag <- diagnFun(NP_multiMod.res)
NP_multiMod.diag #7
NP_multiMod.fail <- failedMdls(NP_multiMod.diag)

for(mdl in NP_multiMod.fail$model){
  if(!is.na(str_extract(mdl, "NP"))){
    mdl1 <- str_remove(mdl, "NP_")
  }
  cat("Re-Running ", mdl)
  multiMdl <- multiMdls[which(multiMdls$name == mdl1), ]
  NP_multiMod.res[mdl] <- multiResponseModels(multiMdl, randForm = "NP_ran.frml")
}

# Save image
save.image(file = "data/RData/plmm_multivariate_NP_results_diagn.RData")
print("====> data/RData/plmm_multivariate_NP_results_diagn.RData")
toc()

remove(list = ls())