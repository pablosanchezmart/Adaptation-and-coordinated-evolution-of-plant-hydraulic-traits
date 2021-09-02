#### MANUSCRIPT TABLES ###################################################################################################################

# P Sanchez-Martínez

print("Preparing manuscript tables...")
remove(list = ls())
source("code/manuscript/1_init.R")

print("loading multivariate models results...")
load("data/RData/plmm_multivariate_NP_results_diagn.RData")
load("data/RData/plmm_multivariate_results_diagn.RData")

print("loading univariate models results...")
load("data/RData/plmm_univariate_NP_results_diagn.RData")
load("data/RData/plmm_univariate_results_diagn.RData")

print("loading PCA results...")
load("data/RData/PCA.RData")

#### FUNCTIONS ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ ####

## Variables number of observations
VarObsFun <- function(ht = ht, vars){
  df <- data.frame("Variable" = vars, "Observations" = NA, "Angiosperms" = NA, "Gymnosperms" = NA)
  for(var in vars){
    df[df$Variable == var, "Observations"] <- length(which(!is.na(ht[, var])))
    df[df$Variable == var, "Angiosperms"] <-  length(which(!is.na(ht[ht$group == "Angiosperms", var])))
    df[df$Variable == var, "Gymnosperms"] <-  length(which(!is.na(ht[ht$group == "Gymnosperms", var])))
  }
  return(df)
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

bestMultiMdlsSumFun <- function(){
  rslts <- list()
  evoCor.res <- data.frame()
  phyInd.res <- data.frame()
  evoCor.df <- read.csv(paste0(output.dir, "EvoCorr_multiModels.csv"), header = T)
  phyInd.df <- read.csv(paste0(output.dir, "phylo_signal_multiModels.csv"), header = T)
  for(i in 1:length(bestMultiMdls$type)){
    vars <- unlist(strsplit(as.character(bestMultiMdls[i, "variables"]), ", "))[c(1:2)]
    vrs <- paste0(vars[1], "_", vars[2], ", PC1, PC2, PC3")
    mdl.name <- as.character(paste0(as.character(bestMultiMdls[i, "type"]), "_", as.character(bestMultiMdls[i, "group"]), "_", vrs))
    for(j in 1:length(multiMod.res)){
    if(multiMod.res[[j]]$name == mdl.name){
      smry <- multiMod.res[[j]]$summary
      smry <- smry$solutions
      NumObs <- multiMod.res[[j]]$NumObs
      rslts[[mdl.name]]$summary <- smry
      rslts[[mdl.name]]$NumObs <- NumObs
      write.csv(smry, paste0(manuscript.dir, "tables/bestMultiModels/", mdl.name, ".csv"), row.names = T)
      print(paste0("====>",manuscript.dir, "tables/bestMultiModels/", mdl.name, ".csv"))
      # PhyloIndexes best models
      phyInd <- unique(phyInd.df[which(phyInd.df$model == mdl.name), ])
      phyInd$model <- str_replace_all(phyInd$model, ", ", "_")
      phyInd$NumObs <- NumObs
      phyInd.res <- rbind(phyInd.res, phyInd)
      # EVoCor best models
      evoCor <- evoCor.df[which(evoCor.df$model == mdl.name), ]
      evoCor$model <- str_replace_all(evoCor$model, ", ", "_")
      evoCor.res <- rbind(evoCor.res, evoCor) 
  } else {
    next()
  }
  }
  }
  rslts$evoCor <- evoCor.res
  rslts$phyInd.res <- phyInd.res
  write.csv(evoCor.res, paste0(manuscript.dir, "tables/bestMultiModels/evoCor.csv"), row.names = F)
  print(paste0("====>",manuscript.dir, "tables/bestMultiModels/evoCor.csv"))
  write.csv(phyInd.res, paste0(manuscript.dir, "tables/bestMultiModels/phyloIndexes.csv"), row.names = F)
  print(paste0("====>",manuscript.dir, "tables/bestMultiModels/phyloIndexes.csv"))
  return(rslts)
}

#### VARIANCE PARTITION (Table 1) --------------------------------------------------------------------------------------------------####

uni.df <- read.csv(paste0(output.dir, "phylo_signal_uniModels.csv"), header = T)
write.csv(uni.df, paste0(manuscript.dir, "tables/variancePartition/variance_partition_uniModels.csv"), row.names = F)
print(paste("==>", manuscript.dir, "tables/variancePartition/variance_partition_uniModels.csv"))

### BEST MULTI-MODELS TABLE (Table 2)---------------------------------------------------------------------------------------------####

bestMultiMdls <- read.csv(paste0(manuscript.dir, "tables/bestMultiModels/bestModelsList.csv"), header = T)
bestMultiMdlsSum <- bestMultiMdlsSumFun()
print(paste("==>", manuscript.dir, "tables/supplementaryTables/variance_partition_NP_uniModels.csv"))

#### VARIABLES OBSERVATIONS (Table S1) ----------------------------------------------------------------------------------------------####

vars <- c(envVars, hydVars, "HSM")
VarObs <- VarObsFun(ht, vars)
write.csv(VarObs, paste0(manuscript.dir, "tables/supplementaryTables/N_Observations.csv"), row.names = F)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/N_Observations.csv"))

Nang <- length(which(!is.na(ht[ht$group == "Angiosperms", "Species"])))
Ngym <- length(which(!is.na(ht[ht$group == "Gymnosperms", "Species"])))

cat("Number of Angiosperms: ", Nang, "\nNumber of gymnosperms: ", Ngym)

#### CONTRIBUTIONS TO THE PCs (Table S2) --------------------------------------------------------------------------------------------####

varContr <- res.var$contrib[, 1:3]
write.csv(varContr, paste0(manuscript.dir, "tables/supplementaryTables/Contributions_to_PCs.csv"), row.names = T)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/Contributions_to_PCs.csv"))

varCor <- res.var$cor[, 1:3]*-1
write.csv(varCor, paste0(manuscript.dir, "tables/supplementaryTables/Correlations_with_PCs.csv"), row.names = T)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/Correlations_with_PCs.csv"))

#### NON-PHYLOGENETIC VARIANCE PARTITION (Table S4)--------------------------------------------------------------------------------####

uni.df <- read.csv(paste0(output.dir, "NP_variances_uniModels.csv"), header = T)
write.csv(uni.df, paste0(manuscript.dir, "tables/supplementaryTables/variance_partition_NP_uniModels.csv"), row.names = F)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/variance_partition_NP_uniModels.csv"))

#### SIGNIFICANT EVOLUTIONARY CORRELATIONS (Table S5) -----------------------------------------------------------------------------------####

# Significant evolutionary correlations
signEvoCor.df <- read.csv(paste0(output.dir, "SignEvoCorr_multiModels.csv"), header = T)
write.csv(signEvoCor.df, paste0(manuscript.dir, "tables/supplementaryTables/Significant_EvolutionaryCorrelations.csv"), row.names = F)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/Significant_EvolutionaryCorrelations.csv"))

#### UNIVARIATE MODELS DESCRIPTION (Table S6) --------------------------------------------------------------------------------------------------####

uni.df <- read.csv(paste0(output.dir, "descriptionAndDIC_uniModels.csv"), header = T)

r2.df <- read.csv(paste0(output.dir, "R2_uniModels.csv"), header = T)
r2.df <- arrange(r2.df, model)
r2_np.df <- read.csv(paste0(output.dir, "NP_R2_uniModels.csv"), header = T)
r2_np.df <- arrange(r2_np.df, model)

r2 <- cbind(r2.df, "NP_R2m_mean" = r2_np.df$R2m_mean, "NP_R2m_lowerHPD" = r2_np.df$R2m_lowerHPD, "NP_R2m_upperHPD" = r2_np.df$R2m_upperHPD,
            "NP_R2c_mean" = r2_np.df$R2c_mean, "NP_R2c_lowerHPD" = r2_np.df$R2c_lowerHPD, "NP_R2c_upper" = r2_np.df$R2c_upperHPD)
r2$name <- r2$model

df <- merge(uni.df, r2, by = c("name"), suffixes = "")
df <- arrange(df, variables, DIC)
df <- df %>% select(name, fix.frml, DIC, R2m_mean, R2c_mean,  NP_R2m_mean, NP_R2c_mean)

write.csv(df, paste0(manuscript.dir, "tables/supplementaryTables/uniModels_description.csv"), row.names = F)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/uniModels_description.csv"))

#### MULTIVARIATE MODELS EVOLUTIONARY CORRELATIONS (Tables S7) -----------------------------------------------------------------------------------####

multi.df <- read.csv(paste0(output.dir, "description_multiModels.csv"), header = T)

# All conbinations and models

evoCor.df <- read.csv(paste0(output.dir, "EvoCorr_multiModels.csv"), header = T)

evoCor.df$DIC <- numeric(length(evoCor.df$model))
evoCor.df$fix.frml <- vector(mode = "character", length = length(evoCor.df$model))
evoCor.df$model <- as.character(evoCor.df$model)

for(i in 1:length(evoCor.df$model)){
  evoCor.df$DIC[i] <- multi.df[which(multi.df$name == evoCor.df$model[i]), "DIC"]
  evoCor.df$fix.frml[i] <- as.character(multi.df[which(multi.df$name == evoCor.df$model[i]), "fix.frml"])
  evoCor.df$fix.frml[i] <- str_replace(evoCor.df$fix.frml[i], ", ", "_")
}
evoCor.df <- arrange(evoCor.df, var1, var2, DIC)
evoCor.df$model <- NULL

write.csv(evoCor.df, paste0(manuscript.dir, "tables/supplementaryTables/EvolutionaryCorrelations.csv"), row.names = F)
print(paste("==>", manuscript.dir, "tables/supplementaryTables/EvolutionaryCorrelations.csv"))


#### SPECIES LEVEL PHYLOGENY RESULTS (APPENDIX 2) -------------------------------------------------------------------------------------------------------------------------------------------------- ####

uni.df <- read.csv(paste0(output.dir, "extension/Spp_phylo_signal_uniModels.csv"), header = T)
write.csv(uni.df, paste0(manuscript.dir, "tables/SA2/Spp_var_partition.csv"), row.names = F)

multi.df <- read.csv(paste0(output.dir, "extension/Spp_EvoCorr_multiModels.csv"), header = T)
write.csv(multi.df, paste0(manuscript.dir, "TABLES/SA2/Spp_EvoCorr_multiModels.csv"),  row.names = F)
          

#### GENUS LEVEL PHYLOGENY SPECIES LEVEL PHYLOGENY DIFFERENCES CHECKING (APPENDIX 3) -------------------------------------------------------------------------------------------------------------------------------------------------- ####

df <- read.csv(paste0(output.dir, "SignEvoCorr_multiModels_sppInt.csv"), header = T)
write.csv(df, paste0(manuscript.dir, "tables/SA3/SignEvoCorr_multiModels_sppInt.csv"), row.names = F)

df <- read.csv(paste0(output.dir, "SignEvoCorr_multiModels_sppTopo.csv"), header = T)
write.csv(df, paste0(manuscript.dir, "tables/SA3/SignEvoCorr_multiModels_sppTopo.csv"), row.names = F)
