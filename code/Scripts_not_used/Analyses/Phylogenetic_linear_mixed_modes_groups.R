################################# PHYLOGENETIC LINEAR MIXED MODELS BY GROUPS ##################################

# Pablo Sanchez
# 02-05-2019

#http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11

### PACKAGES ###
library(ape)
library(MCMCglmm)

setwd("~/Pablo/data")

### FUNCTIONS ####

completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

# Univariate mcmcglmm Function
mcmcglmmFun <- function(dta, fix.frml, ran.frml, prior, nitte = nitter, burn = burni, thinn = thinni){
  # Linear relationship between respnse and cofactor, with a random effect phylo corresponding the phylogenetic effect. The argument ginverse allows us to include a custom matrix for our random effect phylo, 
  mdl <- MCMCglmm(as.formula(fix.frml), random = as.formula(ran.frml),
                  family = "gaussian", pedigree=phylo, prior = prior,
                  data = dta, nitt = nitter, burni = burn, thin = thinni) 
  rslts <- list()
  # Model
  rslts$model <- mdl
  # Autocorrelation, effective size
  rslts$autocFix <- autocorr.diag(mdl$Sol)
  rslts$autocRan <- autocorr.diag(mdl$VCV)
  rslts$heidelFix <- heidel.diag(mdl$Sol)
  rslts$heidelRan <- heidel.diag(mdl$VCV)
  rslts$effSizeFix <- effectiveSize(mdl$VCV)
  rslts$effSizeRan <- effectiveSize(mdl$Sol)
  ## Phylogenetic signal
  lambda <- mdl$VCV[,'animal']/
    (mdl$VCV[,'animal'] + mdl$VCV[,'genus.Rand'] + mdl$VCV[,'units'])
  rslts$lambda <- lambda
  rslts$lambdaMean <- mean(lambda)
  rslts$lambdaPost.mode <- posterior.mode(lambda)
  rslts$lambdaHPDinterval <- HPDinterval(lambda) # 95% credible interval 
  ## inter-genus
  inter <- mdl$VCV[,'genus.Rand']/
    (mdl$VCV[,'animal'] + mdl$VCV[,'genus.Rand'] + mdl$VCV[,'units'])
  rslts$inter <- inter
  rslts$interMean <- mean(inter)
  rslts$interPost.mode <- posterior.mode(inter)
  rslts$interHPDinterval <- HPDinterval(inter) # 95% credible interval 
  ## Intra-genus
  intra <- mdl$VCV[,'units']/
    (mdl$VCV[,'animal'] + mdl$VCV[,'genus.Rand'] + mdl$VCV[,'units'])
  rslts$intra <- intra
  rslts$intraMean <- mean(intra)
  rslts$intraPost.mode <- posterior.mode(intra)
  rslts$intraHPDinterval <- HPDinterval(intra) # 95% credible interval 
  return(rslts)
}

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
  rslts$R2m_HPD <- HPDinterval(R2m)
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
  rslts$R2c_HPD <- HPDinterval(R2c)
  return(rslts)
}

### DATA ####

# Species datset
ht <- read.csv("HydraTRY_Pablo_2019.csv", header = T)  # Includes three principal components

# Genus tree
ht.tree <- read.tree("hydratry_genus.tre") # Genus

# R data
load("RData/not_used/Phylogenetic_linear_mixed_multi_models_groups.RData")

# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% ht$genus])
plmm.data <- ht[ht$genus %in% phylo$tip.label, ]
length(plmm.data$genus)
plmm.data$family <- NULL # gives problems 
length(unique(plmm.data$genus))

# Genus random effect column
plmm.data$genus.Rand <- plmm.data$genus
plmm.data$animal <- plmm.data$genus

### PARAMETERS ####

## Prior distributions for the fixed and random effects
# G for the random effects and R for the residual variance, corresponds to a inverse-Gamma distribution with shape and scape parameter equal to 0.01, which is relatively canonical. 

# Phylogenetic mixed models
prior2 <- list(
  R = list(V = 1, nu = 0.002), 
  G = list(G1 = list(V = 1, nu = 0.002), 
           G2 = list(V = 1, nu = 0.002))) # two random effects, phylogeny and genus

# With parameter expansion for random effects
priorEx <- list(
  R=list(V = 1, nu = 0.002), 
  G=list(G1=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1),
         G2 = list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1)))
 
 # Cauchy
 priorCau <- list(
   R = list(V = 1, nu = 0.002),
   G = list(G1 = list(V= 1, nu = 1, alpha.mu = 0),
            G2 = (G1 = list(V= 1, nu = 1, alpha.mu = 0))))

   
# Phylogenetic mixed multi models
priorMulti <- list(
  R=list(V=diag(2)/2,nu=2), 
  G=list(G1=list(V=diag(2)/2,nu=2), 
         G2 = list(V=diag(2)/2,nu=2)))

priorMulti2 <- list(
  R=list(V=diag(2)*(0.002/1.002),nu=1.002), 
  G=list(G1=list(V=diag(2)*(0.002/1.002),nu=1.002),
         G2 = list(V=diag(2)/2,nu=2)))

priorMultiCVstr <- list(
  R=list(V=diag(3)*(0.002/1.002),nu=1.002), 
  G=list(G1=list(V=diag(3)*(0.002/1.002),nu=1.002),
         G2 = list(V=diag(3)/2,nu=2)))
# Model parameters
nitter = 8000000
burni = 1000000
thinni = 4000
set.seed(123)

########### PHYLOGENETIC SIGNAL #####################################################################################################

### PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "genus"))
groups <- c("Angiosperms", "Gymnosperms")
ps.groups <- data.frame()
inter.groups <- data.frame()
intra.groups <- data.frame()

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  PC1.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "PC1 ~ 1", ran.frml = "~ animal + genus.Rand", 
                               prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(PC1.groupsMod$model)
  PC1.groupsMod$autocFix
  PC1.groupsMod$autocRan
  PC1.groupsMod$heidelFix
  PC1.groupsMod$heidelRan
  PC1.groupsMod$effSizeFix
  PC1.groupsMod$effSizeRan
  
  # Lambda
  # plot(PC1.groupsMod$lambda)
  PC1.groupsMod$lambdaMean
  PC1.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(PC1.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "PC1"
  ps.group[gr, "Mean"] <- PC1.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- PC1.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- PC1.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "PC1"
  inter.group[gr, "Mean"] <- PC1.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- PC1.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- PC1.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "PC1"
  intra.group[gr, "Mean"] <- PC1.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- PC1.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- PC1.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}


### PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  PC2.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "PC2 ~ 1", ran.frml = "~ animal + genus.Rand", 
                               prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(PC2.groupsMod$model)
  PC2.groupsMod$autocFix
  PC2.groupsMod$autocRan
  PC2.groupsMod$heidelFix
  PC2.groupsMod$heidelRan
  PC2.groupsMod$effSizeFix
  PC2.groupsMod$effSizeRan
  
  # Lambda
  # plot(PC2.groupsMod$lambda)
  PC2.groupsMod$lambdaMean
  PC2.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(PC2.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "PC2"
  ps.group[gr, "Mean"] <- PC2.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- PC2.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- PC2.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "PC2"
  inter.group[gr, "Mean"] <- PC2.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- PC2.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- PC2.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "PC2"
  intra.group[gr, "Mean"] <- PC2.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- PC2.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- PC2.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}


### PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  PC3.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "PC3 ~ 1", ran.frml = "~ animal + genus.Rand", 
                               prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(PC3.groupsMod$model)
  PC3.groupsMod$autocFix
  PC3.groupsMod$autocRan
  PC3.groupsMod$heidelFix
  PC3.groupsMod$heidelRan
  PC3.groupsMod$effSizeFix
  PC3.groupsMod$effSizeRan
  
  # Lambda
  # plot(PC3.groupsMod$lambda)
  PC3.groupsMod$lambdaMean
  PC3.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(PC3.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "PC3"
  ps.group[gr, "Mean"] <- PC3.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- PC3.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- PC3.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "PC3"
  inter.group[gr, "Mean"] <- PC3.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- PC3.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- PC3.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "PC3"
  intra.group[gr, "Mean"] <- PC3.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- PC3.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- PC3.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}


### log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_Ks.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_Ks ~ 1", ran.frml = "~ animal + genus.Rand", 
                                  prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_Ks.groupsMod$model)
  log_Ks.groupsMod$autocFix
  log_Ks.groupsMod$autocRan
  log_Ks.groupsMod$heidelFix
  log_Ks.groupsMod$heidelRan
  log_Ks.groupsMod$effSizeFix
  log_Ks.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_Ks.groupsMod$lambda)
  log_Ks.groupsMod$lambdaMean
  log_Ks.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_Ks.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_Ks"
  ps.group[gr, "Mean"] <- log_Ks.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_Ks.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_Ks.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_Ks"
  inter.group[gr, "Mean"] <- log_Ks.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_Ks.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_Ks.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_Ks"
  intra.group[gr, "Mean"] <- log_Ks.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_Ks.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_Ks.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}


### log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_negP50.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_negP50 ~ 1", ran.frml = "~ animal + genus.Rand", 
                                      prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_negP50.groupsMod$model)
  log_negP50.groupsMod$autocFix
  log_negP50.groupsMod$autocRan
  log_negP50.groupsMod$heidelFix
  log_negP50.groupsMod$heidelRan
  log_negP50.groupsMod$effSizeFix
  log_negP50.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_negP50.groupsMod$lambda)
  log_negP50.groupsMod$lambdaMean
  log_negP50.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_negP50.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_negP50"
  ps.group[gr, "Mean"] <- log_negP50.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_negP50.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_negP50.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_negP50"
  inter.group[gr, "Mean"] <- log_negP50.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_negP50.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_negP50.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_negP50"
  intra.group[gr, "Mean"] <- log_negP50.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_negP50.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_negP50.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}

### log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negMinWP_md", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_negMinWP_md.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_negMinWP_md ~ 1", ran.frml = "~ animal + genus.Rand", 
                                           prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_negMinWP_md.groupsMod$model)
  log_negMinWP_md.groupsMod$autocFix
  log_negMinWP_md.groupsMod$autocRan
  log_negMinWP_md.groupsMod$heidelFix
  log_negMinWP_md.groupsMod$heidelRan
  log_negMinWP_md.groupsMod$effSizeFix
  log_negMinWP_md.groupsMod$effSizeRan
  
  # Lambda
  plot(log_negMinWP_md.groupsMod$lambda)
  log_negMinWP_md.groupsMod$lambdaMean
  log_negMinWP_md.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_negMinWP_md.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_negMinWP_md"
  ps.group[gr, "Mean"] <- log_negMinWP_md.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_negMinWP_md.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_negMinWP_md.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_negMinWP_md"
  inter.group[gr, "Mean"] <- log_negMinWP_md.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_negMinWP_md.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_negMinWP_md.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_negMinWP_md"
  intra.group[gr, "Mean"] <- log_negMinWP_md.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_negMinWP_md.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_negMinWP_md.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}

### log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_Hv.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_Hv ~ 1", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_Hv.groupsMod$model)
  log_Hv.groupsMod$autocFix
  log_Hv.groupsMod$autocRan
  log_Hv.groupsMod$heidelFix
  log_Hv.groupsMod$heidelRan
  log_Hv.groupsMod$effSizeFix
  log_Hv.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_Hv.groupsMod$lambda)
  log_Hv.groupsMod$lambdaMean
  log_Hv.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_Hv.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_Hv"
  ps.group[gr, "Mean"] <- log_Hv.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_Hv.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_Hv.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_Hv"
  inter.group[gr, "Mean"] <- log_Hv.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_Hv.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_Hv.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_Hv"
  intra.group[gr, "Mean"] <- log_Hv.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_Hv.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_Hv.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}

### log_negPItlp ###

# Complete data
data <- completeFun(plmm.data, c("log_negPItlp", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_negPItlp.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_negPItlp ~ 1", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_negPItlp.groupsMod$model)
  log_negPItlp.groupsMod$autocFix
  log_negPItlp.groupsMod$autocRan
  log_negPItlp.groupsMod$heidelFix
  log_negPItlp.groupsMod$heidelRan
  log_negPItlp.groupsMod$effSizeFix
  log_negPItlp.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_negPItlp.groupsMod$lambda)
  log_negPItlp.groupsMod$lambdaMean
  log_negPItlp.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_negPItlp.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_negPItlp"
  ps.group[gr, "Mean"] <- log_negPItlp.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_negPItlp.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_negPItlp.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_negPItlp"
  inter.group[gr, "Mean"] <- log_negPItlp.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_negPItlp.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_negPItlp.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_negPItlp"
  intra.group[gr, "Mean"] <- log_negPItlp.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_negPItlp.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_negPItlp.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}
### log_Kl ###

# Complete data
data <- completeFun(plmm.data, c("log_Kl", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_Kl.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_Kl ~ 1", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_Kl.groupsMod$model)
  log_Kl.groupsMod$autocFix
  log_Kl.groupsMod$autocRan
  log_Kl.groupsMod$heidelFix
  log_Kl.groupsMod$heidelRan
  log_Kl.groupsMod$effSizeFix
  log_Kl.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_Kl.groupsMod$lambda)
  log_Kl.groupsMod$lambdaMean
  log_Kl.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_Kl.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_Kl"
  ps.group[gr, "Mean"] <- log_Kl.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_Kl.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_Kl.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_Kl"
  inter.group[gr, "Mean"] <- log_Kl.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_Kl.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_Kl.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_Kl"
  intra.group[gr, "Mean"] <- log_Kl.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_Kl.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_Kl.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}
### log_Hmax ###

# Complete data
data <- completeFun(plmm.data, c("log_Hmax", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_Hmax.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_Hmax ~ 1", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_Hmax.groupsMod$model)
  log_Hmax.groupsMod$autocFix
  log_Hmax.groupsMod$autocRan
  log_Hmax.groupsMod$heidelFix
  log_Hmax.groupsMod$heidelRan
  log_Hmax.groupsMod$effSizeFix
  log_Hmax.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_Hmax.groupsMod$lambda)
  log_Hmax.groupsMod$lambdaMean
  log_Hmax.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_Hmax.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_Hmax"
  ps.group[gr, "Mean"] <- log_Hmax.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_Hmax.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_Hmax.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_Hmax"
  inter.group[gr, "Mean"] <- log_Hmax.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_Hmax.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_Hmax.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_Hmax"
  intra.group[gr, "Mean"] <- log_Hmax.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_Hmax.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_Hmax.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}
### Rd ###

# Complete data
data <- completeFun(plmm.data, c("Rd", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  Rd.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "Rd ~ 1", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(Rd.groupsMod$model)
  Rd.groupsMod$autocFix
  Rd.groupsMod$autocRan
  Rd.groupsMod$heidelFix
  Rd.groupsMod$heidelRan
  Rd.groupsMod$effSizeFix
  Rd.groupsMod$effSizeRan
  
  # Lambda
  # plot(Rd.groupsMod$lambda)
  Rd.groupsMod$lambdaMean
  Rd.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(Rd.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "Rd"
  ps.group[gr, "Mean"] <- Rd.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- Rd.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- Rd.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "Rd"
  inter.group[gr, "Mean"] <- Rd.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- Rd.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- Rd.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "Rd"
  intra.group[gr, "Mean"] <- Rd.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- Rd.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- Rd.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}
### log_gmax ###

# Complete data
data <- completeFun(plmm.data, c("log_gmax", "genus"))

for(gr in groups){
  data.gr <- data[which(data$group == gr), ]
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data.gr$genus])
  data.gr <- data.gr[data.gr$genus %in% phylo$tip.label, ]
  length(data.gr$genus)
  
  # Model
  log_gmax.groupsMod <- mcmcglmmFun(dta = data.gr, fix.frml = "log_gmax ~ 1", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2, nitter, burni, thinni)
  
  ## Autocorrelation and convergence
  # plot(log_gmax.groupsMod$model)
  log_gmax.groupsMod$autocFix
  log_gmax.groupsMod$autocRan
  log_gmax.groupsMod$heidelFix
  log_gmax.groupsMod$heidelRan
  log_gmax.groupsMod$effSizeFix
  log_gmax.groupsMod$effSizeRan
  
  # Lambda
  # plot(log_gmax.groupsMod$lambda)
  log_gmax.groupsMod$lambdaMean
  log_gmax.groupsMod$lambdaHPDinterval
  
  # Model results
  summary(log_gmax.groupsMod$model)
  
  # Save
  ps.group <- data.frame()
  ps.group[gr, "Variable"] <- "log_gmax"
  ps.group[gr, "Mean"] <- log_gmax.groupsMod$lambdaMean
  ps.group[gr, "lower HDP"] <- log_gmax.groupsMod$lambdaHPDinterval[1]
  ps.group[gr, "higher HDP"] <- log_gmax.groupsMod$lambdaHPDinterval[2]
  ps.groups <- rbind(ps.groups, ps.group)
  
  inter.group <- data.frame()
  inter.group[gr, "Variable"] <- "log_gmax"
  inter.group[gr, "Mean"] <- log_gmax.groupsMod$interMean
  inter.group[gr, "lower HDP"] <- log_gmax.groupsMod$interHPDinterval[1]
  inter.group[gr, "higher HDP"] <- log_gmax.groupsMod$interHPDinterval[2]
  inter.groups <- rbind(inter.groups, inter.group)
  
  intra.group <- data.frame()
  intra.group[gr, "Variable"] <- "log_gmax"
  intra.group[gr, "Mean"] <- log_gmax.groupsMod$intraMean
  intra.group[gr, "lower HDP"] <- log_gmax.groupsMod$intraHPDinterval[1]
  intra.group[gr, "higher HDP"] <- log_gmax.groupsMod$intraHPDinterval[2]
  intra.groups <- rbind(intra.groups, intra.group)
}
ps.groups
inter.groups
intra.groups

########### HT ~ ENVIRONMENT UNIVARIATE #####################################################################################################

### log_Ks ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_PC1_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC1 * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_PC1_gr.mod$model)
log_Ks_PC1_gr.mod$autocFix
log_Ks_PC1_gr.mod$autocRan
log_Ks_PC1_gr.mod$heidelFix
log_Ks_PC1_gr.mod$heidelRan
log_Ks_PC1_gr.mod$effSizeFix
log_Ks_PC1_gr.mod$effSizeRan

# Lambda
plot(log_Ks_PC1_gr.mod$lambda)
log_Ks_PC1_gr.mod$lambdaMean
log_Ks_PC1_gr.mod$lambdaHPDinterval

# Model results
summary(log_Ks_PC1_gr.mod$model)


### log_Ks ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_PC2_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC2 * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_PC2_gr.mod$model)
log_Ks_PC2_gr.mod$autocFix
log_Ks_PC2_gr.mod$autocRan
log_Ks_PC2_gr.mod$heidelFix
log_Ks_PC2_gr.mod$heidelRan
log_Ks_PC2_gr.mod$effSizeFix
log_Ks_PC2_gr.mod$effSizeRan

# Lambda
plot(log_Ks_PC2_gr.mod$lambda)
log_Ks_PC2_gr.mod$lambdaMean
log_Ks_PC2_gr.mod$lambdaHPDinterval

# Model results
summary(log_Ks_PC2_gr.mod$model)

### log_Ks ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_PC3_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ PC3 * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_PC3_gr.mod$model)
log_Ks_PC3_gr.mod$autocFix
log_Ks_PC3_gr.mod$autocRan
log_Ks_PC3_gr.mod$heidelFix
log_Ks_PC3_gr.mod$heidelRan
log_Ks_PC3_gr.mod$effSizeFix
log_Ks_PC3_gr.mod$effSizeRan

# Lambda
plot(log_Ks_PC3_gr.mod$lambda)
log_Ks_PC3_gr.mod$lambdaMean
log_Ks_PC3_gr.mod$lambdaHPDinterval

# Model results
summary(log_Ks_PC3_gr.mod$model)


### log_negP50 ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_PC1_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC1 * group", ran.frml = "~ animal + genus.Rand", 
                                     prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_PC1_gr.mod$model)
log_negP50_PC1_gr.mod$autocFix
log_negP50_PC1_gr.mod$autocRan
log_negP50_PC1_gr.mod$heidelFix
log_negP50_PC1_gr.mod$heidelRan
log_negP50_PC1_gr.mod$effSizeFix
log_negP50_PC1_gr.mod$effSizeRan

# Lambda
plot(log_negP50_PC1_gr.mod$lambda)
log_negP50_PC1_gr.mod$lambdaMean
log_negP50_PC1_gr.mod$lambdaHPDinterval

# Model results
summary(log_negP50_PC1_gr.mod$model)


### log_negP50 ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_PC2_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC2 * group", ran.frml = "~ animal + genus.Rand", 
                                     prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_PC2_gr.mod$model)
log_negP50_PC2_gr.mod$autocFix
log_negP50_PC2_gr.mod$autocRan
log_negP50_PC2_gr.mod$heidelFix
log_negP50_PC2_gr.mod$heidelRan
log_negP50_PC2_gr.mod$effSizeFix
log_negP50_PC2_gr.mod$effSizeRan

# Lambda
plot(log_negP50_PC2_gr.mod$lambda)
log_negP50_PC2_gr.mod$lambdaMean
log_negP50_PC2_gr.mod$lambdaHPDinterval

# Model results
summary(log_negP50_PC2_gr.mod$model)

r2m <- R2mFun(log_negP50_PC2_gr.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negP50_PC2_gr.mod$model)
r2c$R2c_mean
r2c$R2c_HPD


### log_negP50 ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_PC3_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ PC3 * group", ran.frml = "~ animal + genus.Rand", 
                                     prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_PC3_gr.mod$model)
log_negP50_PC3_gr.mod$autocFix
log_negP50_PC3_gr.mod$autocRan
log_negP50_PC3_gr.mod$heidelFix
log_negP50_PC3_gr.mod$heidelRan
log_negP50_PC3_gr.mod$effSizeFix
log_negP50_PC3_gr.mod$effSizeRan

# Lambda
plot(log_negP50_PC3_gr.mod$lambda)
log_negP50_PC3_gr.mod$lambdaMean
log_negP50_PC3_gr.mod$lambdaHPDinterval

# Model results
summary(log_negP50_PC3_gr.mod$model)


### log_negMinWP_md ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_PC1_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC1 * group", ran.frml = "~ animal + genus.Rand", 
                                          prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_PC1_gr.mod$model)
log_negMinWP_md_PC1_gr.mod$autocFix
log_negMinWP_md_PC1_gr.mod$autocRan
log_negMinWP_md_PC1_gr.mod$heidelFix
log_negMinWP_md_PC1_gr.mod$heidelRan
log_negMinWP_md_PC1_gr.mod$effSizeFix
log_negMinWP_md_PC1_gr.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_PC1_gr.mod$lambda)
log_negMinWP_md_PC1_gr.mod$lambdaMean
log_negMinWP_md_PC1_gr.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_PC1_gr.mod$model)

r2m <- R2mFun(log_negMinWP_md_PC1_gr.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negMinWP_md_PC1_gr.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negMinWP_md ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_PC2_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC2 * group", ran.frml = "~ animal + genus.Rand", 
                                          prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_PC2_gr.mod$model)
log_negMinWP_md_PC2_gr.mod$autocFix
log_negMinWP_md_PC2_gr.mod$autocRan
log_negMinWP_md_PC2_gr.mod$heidelFix
log_negMinWP_md_PC2_gr.mod$heidelRan
log_negMinWP_md_PC2_gr.mod$effSizeFix
log_negMinWP_md_PC2_gr.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_PC2_gr.mod$lambda)
log_negMinWP_md_PC2_gr.mod$lambdaMean
log_negMinWP_md_PC2_gr.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_PC2_gr.mod$model)

### log_negMinWP_md ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_PC3_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ PC3 * group", ran.frml = "~ animal + genus.Rand", 
                                          prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_PC3_gr.mod$model)
log_negMinWP_md_PC3_gr.mod$autocFix
log_negMinWP_md_PC3_gr.mod$autocRan
log_negMinWP_md_PC3_gr.mod$heidelFix
log_negMinWP_md_PC3_gr.mod$heidelRan
log_negMinWP_md_PC3_gr.mod$effSizeFix
log_negMinWP_md_PC3_gr.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_PC3_gr.mod$lambda)
log_negMinWP_md_PC3_gr.mod$lambdaMean
log_negMinWP_md_PC3_gr.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_PC3_gr.mod$model)

### log_Hv ~ PC1 ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_PC1_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC1 * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_PC1_gr.mod$model)
log_Hv_PC1_gr.mod$autocFix
log_Hv_PC1_gr.mod$autocRan
log_Hv_PC1_gr.mod$heidelFix
log_Hv_PC1_gr.mod$heidelRan
log_Hv_PC1_gr.mod$effSizeFix
log_Hv_PC1_gr.mod$effSizeRan

# Lambda
plot(log_Hv_PC1_gr.mod$lambda)
log_Hv_PC1_gr.mod$lambdaMean
log_Hv_PC1_gr.mod$lambdaHPDinterval

# Model results
summary(log_Hv_PC1_gr.mod$model)

### log_Hv ~ PC2 ###

# Complete data
data <- completeFun(plmm.data, c("PC2", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_PC2_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC2 * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_PC2_gr.mod$model)
log_Hv_PC2_gr.mod$autocFix
log_Hv_PC2_gr.mod$autocRan
log_Hv_PC2_gr.mod$heidelFix
log_Hv_PC2_gr.mod$heidelRan
log_Hv_PC2_gr.mod$effSizeFix
log_Hv_PC2_gr.mod$effSizeRan

# Lambda
plot(log_Hv_PC2_gr.mod$lambda)
log_Hv_PC2_gr.mod$lambdaMean
log_Hv_PC2_gr.mod$lambdaHPDinterval

# Model results
summary(log_Hv_PC2_gr.mod$model)


### log_Hv ~ PC3 ###

# Complete data
data <- completeFun(plmm.data, c("PC3", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_PC3_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ PC3 * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_PC3_gr.mod$model)
log_Hv_PC3_gr.mod$autocFix
log_Hv_PC3_gr.mod$autocRan
log_Hv_PC3_gr.mod$heidelFix
log_Hv_PC3_gr.mod$heidelRan
log_Hv_PC3_gr.mod$effSizeFix
log_Hv_PC3_gr.mod$effSizeRan

# Lambda
plot(log_Hv_PC3_gr.mod$lambda)
log_Hv_PC3_gr.mod$lambdaMean
log_Hv_PC3_gr.mod$lambdaHPDinterval

# Model results
summary(log_Hv_PC3_gr.mod$model)


########### HT ~ ENVIRONMENT #####################################################################################################

### log_Ks ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Ks"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_env_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_env_gr.mod$model)
log_Ks_env_gr.mod$autocFix
log_Ks_env_gr.mod$autocRan
log_Ks_env_gr.mod$heidelFix
log_Ks_env_gr.mod$heidelRan
log_Ks_env_gr.mod$effSizeFix
log_Ks_env_gr.mod$effSizeRan

# Lambda
plot(log_Ks_env_gr.mod$lambda)
log_Ks_env_gr.mod$lambdaMean
log_Ks_env_gr.mod$lambdaHPDinterval

# Model results
summary(log_Ks_env_gr.mod$model)
R2cFun(log_Ks_env_gr.mod$model)

### log_negP50 ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negP50"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negP50_env_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ animal + genus.Rand", 
                                     prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_env_gr.mod$model)
log_negP50_env_gr.mod$autocFix
log_negP50_env_gr.mod$autocRan
log_negP50_env_gr.mod$heidelFix
log_negP50_env_gr.mod$heidelRan
log_negP50_env_gr.mod$effSizeFix
log_negP50_env_gr.mod$effSizeRan

# Lambda
plot(log_negP50_env_gr.mod$lambda)
log_negP50_env_gr.mod$lambdaMean
log_negP50_env_gr.mod$lambdaHPDinterval

# Model results
summary(log_negP50_env_gr.mod$model)
R2cFun(log_negP50_env_gr.mod$model)

### log_negMinWP_md ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_negMinWP_md"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_negMinWP_md_env_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ animal + genus.Rand", 
                                          prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_env_gr.mod$model)
log_negMinWP_md_env_gr.mod$autocFix
log_negMinWP_md_env_gr.mod$autocRan
log_negMinWP_md_env_gr.mod$heidelFix
log_negMinWP_md_env_gr.mod$heidelRan
log_negMinWP_md_env_gr.mod$effSizeFix
log_negMinWP_md_env_gr.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_env_gr.mod$lambda)
log_negMinWP_md_env_gr.mod$lambdaMean
log_negMinWP_md_env_gr.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_env_gr.mod$model)
R2cFun(log_negMinWP_md_env_gr.mod$model)

### log_Hv ~ environment ###

# Complete data
data <- completeFun(plmm.data, c("PC1", "log_Hv"))
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_env_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ (PC1 + PC2 + PC3) * group", ran.frml = "~ animal + genus.Rand", 
                                 prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_env_gr.mod$model)
log_Hv_env_gr.mod$autocFix
log_Hv_env_gr.mod$autocRan
log_Hv_env_gr.mod$heidelFix
log_Hv_env_gr.mod$heidelRan
log_Hv_env_gr.mod$effSizeFix
log_Hv_env_gr.mod$effSizeRan

# Lambda
plot(log_Hv_env_gr.mod$lambda)
log_Hv_env_gr.mod$lambdaMean
log_Hv_env_gr.mod$lambdaHPDinterval

# Model results
summary(log_Hv_env_gr.mod$model)
R2cFun(log_Hv_env_gr.mod$model)

########### HT ~ HT ##############################################################################################################


### log_Ks ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ log_negP50 * group", ran.frml = "~ animal + genus.Rand",   # us(group):animal + us(group):genus.Rand (allow to include the interaction)
                                        prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_log_negP50_gr.mod$model)
log_Ks_log_negP50_gr.mod$autocFix
log_Ks_log_negP50_gr.mod$autocRan
log_Ks_log_negP50_gr.mod$heidelFix
log_Ks_log_negP50_gr.mod$heidelRan
log_Ks_log_negP50_gr.mod$effSizeFix
log_Ks_log_negP50_gr.mod$effSizeRan

# Lambda
plot(log_Ks_log_negP50_gr.mod$lambda)
log_Ks_log_negP50_gr.mod$lambdaMean
log_Ks_log_negP50_gr.mod$lambdaHPDinterval

# Model results
summary(log_Ks_log_negP50_gr.mod$model)


### log_negP50 ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Ks", "log_negP50"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negP50_log_Ks_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ log_Ks * group", ran.frml = "~ animal + genus.Rand",   # us(group):animal + us(group):genus.Rand (allow to include the interaction)
                                        prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_log_Ks_gr.mod$model)
log_negP50_log_Ks_gr.mod$autocFix
log_negP50_log_Ks_gr.mod$autocRan
log_negP50_log_Ks_gr.mod$heidelFix
log_negP50_log_Ks_gr.mod$heidelRan
log_negP50_log_Ks_gr.mod$effSizeFix
log_negP50_log_Ks_gr.mod$effSizeRan

# Lambda
plot(log_negP50_log_Ks_gr.mod$lambda)
log_negP50_log_Ks_gr.mod$lambdaMean
log_negP50_log_Ks_gr.mod$lambdaHPDinterval

# Model results
summary(log_negP50_log_Ks_gr.mod$model)


### log_negP50 ~ log_negMinWP_md ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negP50_log_negMinWP_md_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negP50 ~ log_negMinWP_md * group", ran.frml = "~ animal + genus.Rand", 
                                                 prior = prior2)

## Autocorrelation and convergence
plot(log_negP50_log_negMinWP_md_gr.mod$model)
log_negP50_log_negMinWP_md_gr.mod$autocFix
log_negP50_log_negMinWP_md_gr.mod$autocRan
log_negP50_log_negMinWP_md_gr.mod$heidelFix
log_negP50_log_negMinWP_md_gr.mod$heidelRan
log_negP50_log_negMinWP_md_gr.mod$effSizeFix
log_negP50_log_negMinWP_md_gr.mod$effSizeRan

# Lambda
plot(log_negP50_log_negMinWP_md_gr.mod$lambda)
log_negP50_log_negMinWP_md_gr.mod$lambdaMean
log_negP50_log_negMinWP_md_gr.mod$lambdaHPDinterval

# Model results
summary(log_negP50_log_negMinWP_md_gr.mod$model)

r2m <- R2mFun(log_negP50_log_negMinWP_md_gr.mod$model)
r2m$R2m_mean
r2m$R2m_HPD
r2c <- R2cFun(log_negP50_log_negMinWP_md_gr.mod$model)
r2c$R2c_mean
r2c$R2c_HPD

### log_negMinWP_md ~ log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
#Model
log_negMinWP_md_log_negP50_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_negMinWP_md ~ log_negP50 * group", ran.frml = "~ animal + genus.Rand", 
                                                 prior = prior2)

## Autocorrelation and convergence
plot(log_negMinWP_md_log_negP50_gr.mod$model)
log_negMinWP_md_log_negP50_gr.mod$autocFix
log_negMinWP_md_log_negP50_gr.mod$autocRan
log_negMinWP_md_log_negP50_gr.mod$heidelFix
log_negMinWP_md_log_negP50_gr.mod$heidelRan
log_negMinWP_md_log_negP50_gr.mod$effSizeFix
log_negMinWP_md_log_negP50_gr.mod$effSizeRan

# Lambda
plot(log_negMinWP_md_log_negP50_gr.mod$lambda)
log_negMinWP_md_log_negP50_gr.mod$lambdaMean
log_negMinWP_md_log_negP50_gr.mod$lambdaHPDinterval

# Model results
summary(log_negMinWP_md_log_negP50_gr.mod$model)


### log_Hv ~ log_Ks ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Hv_log_Ks_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Hv ~ log_Ks * group", ran.frml = "~ animal + genus.Rand", 
                                    prior = prior2)

## Autocorrelation and convergence
plot(log_Hv_log_Ks_gr.mod$model)
log_Hv_log_Ks_gr.mod$autocFix
log_Hv_log_Ks_gr.mod$autocRan
log_Hv_log_Ks_gr.mod$heidelFix
log_Hv_log_Ks_gr.mod$heidelRan
log_Hv_log_Ks_gr.mod$effSizeFix
log_Hv_log_Ks_gr.mod$effSizeRan

# Lambda
plot(log_Hv_log_Ks_gr.mod$lambda)
log_Hv_log_Ks_gr.mod$lambdaMean
log_Hv_log_Ks_gr.mod$lambdaHPDinterval

# Model results
summary(log_Hv_log_Ks_gr.mod$model)


### log_Ks ~ log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))  # more than 100 NAs, this could be 
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)
# Model
log_Ks_log_Hv_gr.mod <- mcmcglmmFun(dta = data, fix.frml = "log_Ks ~ log_Hv * group", ran.frml = "~ animal + genus.Rand", 
                                    prior = prior2)

## Autocorrelation and convergence
plot(log_Ks_log_Hv_gr.mod$model)
log_Ks_log_Hv_gr.mod$autocFix
log_Ks_log_Hv_gr.mod$autocRan
log_Ks_log_Hv_gr.mod$heidelFix
log_Ks_log_Hv_gr.mod$heidelRan
log_Ks_log_Hv_gr.mod$effSizeFix
log_Ks_log_Hv_gr.mod$effSizeRan

# Lambda
plot(log_Ks_log_Hv_gr.mod$lambda)
log_Ks_log_Hv_gr.mod$lambdaMean
log_Ks_log_Hv_gr.mod$lambdaHPDinterval

# Model results
summary(log_Ks_log_Hv_gr.mod$model)


############################# HT + HT ~ -1 + group ##############################################################################

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
meanKsP50_gr.multiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:group, 
                                  random = ~us(trait):animal + us(trait):genus.Rand, 
                                  rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                  prior = priorMulti2, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsP50_gr.multiMod$Sol)
plot(meanKsP50_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsP50_gr.multiMod$Sol)
autocorr.diag(meanKsP50_gr.multiMod$VCV)
# Convergence
heidel.diag(meanKsP50_gr.multiMod$Sol)
heidel.diag(meanKsP50_gr.multiMod$VCV)

## Phylogenetic signal

# log_Ks
physig_log_Ks_gr.multiMod <- meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]+ meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"])
mean(physig_log_Ks_gr.multiMod)
HPDinterval(physig_log_Ks_gr.multiMod, 0.95)

# log_negP50
physig_log_negP50_gr.multiMod <- meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]+ meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"])
mean(physig_log_negP50_gr.multiMod)
HPDinterval(physig_log_negP50_gr.multiMod, 0.95)

## Inter-genus

# log_Ks
inter_log_Ks_gr.multiMod <- meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]+ meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"])
mean(inter_log_Ks_gr.multiMod)
HPDinterval(inter_log_Ks_gr.multiMod, 0.95)

# log_negP50
inter_log_negP50_gr.multiMod <- meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]+ meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"])
mean(inter_log_negP50_gr.multiMod)
HPDinterval(inter_log_negP50_gr.multiMod, 0.95)

## Intra-genus

# log_Ks
intra_log_Ks_gr.multiMod <- meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]+ meanKsP50_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"])
mean(intra_log_Ks_gr.multiMod)
HPDinterval(intra_log_Ks_gr.multiMod, 0.95)

# log_negP50
intra_log_negP50_gr.multiMod <- meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"] + meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"])
mean(intra_log_negP50_gr.multiMod)
HPDinterval(intra_log_negP50_gr.multiMod, 0.95)

## Evolutive correlation
cor_log_Ks_log_negP50_gr.multiMod <- meanKsP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_Ks.animal"]/
  sqrt(meanKsP50_gr.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*meanKsP50_gr.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(cor_log_Ks_log_negP50_gr.multiMod)
HPDinterval(cor_log_Ks_log_negP50_gr.multiMod, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
meanKsHv_gr.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1  + trait:group, 
                                 random = ~us(trait):animal + us(trait):genus.Rand, 
                                 rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                 prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanKsHv_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanKsHv_gr.multiMod$Sol)
autocorr.diag(meanKsHv_gr.multiMod$VCV)
# Convergence
heidel.diag(meanKsHv_gr.multiMod$Sol)
heidel.diag(meanKsHv_gr.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks_gr.multiMod <- meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]+ meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"])
mean(physig_log_Ks_gr.multiMod)
HPDinterval(physig_log_Ks_gr.multiMod, 0.95)

# log_Hv
physig_log_Hv_gr.multiMod <- meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"] + meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"])
mean(physig_log_Hv_gr.multiMod)
HPDinterval(physig_log_Hv_gr.multiMod, 0.95)

# Inter-genus

# log_Ks
inter_log_Ks_gr.multiMod <- meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]+ meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"])
mean(inter_log_Ks_gr.multiMod)
HPDinterval(inter_log_Ks_gr.multiMod, 0.95)

# log_Hv
inter_log_Hv_gr.multiMod <- meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"]/ 
  (meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"] + meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"])
mean(inter_log_Hv_gr.multiMod)
HPDinterval(inter_log_Hv_gr.multiMod, 0.95)

# Intra-genus

# log_Ks
intra_log_Ks_gr.multiMod <- meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"] + meanKsHv_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"])
mean(intra_log_Ks_gr.multiMod)
HPDinterval(intra_log_Ks_gr.multiMod, 0.95)

# log_Hv
intra_log_Hv_gr.multiMod <- meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"]/ 
  (meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"] + meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"])
mean(intra_log_Hv_gr.multiMod)
HPDinterval(intra_log_Hv_gr.multiMod, 0.95)

# Evolutive correlation
cor_log_Ks_log_Hv_gr.multiMod <- meanKsHv_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(meanKsHv_gr.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*meanKsHv_gr.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv_gr.multiMod)
HPDinterval(cor_log_Ks_log_Hv_gr.multiMod, 0.95)


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
meanMinWPP50_gr.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1  + trait:group, 
                                     random = ~us(trait):animal + us(trait):genus.Rand, 
                                     rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                     prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(meanMinWPP50_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(meanMinWPP50_gr.multiMod$Sol)
autocorr.diag(meanMinWPP50_gr.multiMod$VCV)
# Convergence
heidel.diag(meanMinWPP50_gr.multiMod$Sol)
heidel.diag(meanMinWPP50_gr.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]+ meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"])
mean(physig_log_negMinWP_md_gr.multiMod)
HPDinterval(physig_log_negMinWP_md_gr.multiMod, 0.95)

# log_negP50
physig_log_negP50_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"] + meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"])
mean(physig_log_negP50_gr.multiMod)
HPDinterval(physig_log_negP50_gr.multiMod, 0.95)

# Inter-genus

# log_negMinWP_md
inter_log_negMinWP_md_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"]/ 
  (meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]+ meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"])
mean(inter_log_negMinWP_md_gr.multiMod)
HPDinterval(inter_log_negMinWP_md_gr.multiMod, 0.95)

# log_negP50
inter_log_negP50_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"] + meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"])
mean(inter_log_negP50_gr.multiMod)
HPDinterval(inter_log_negP50_gr.multiMod, 0.95)

# Intra-genus

# log_negMinWP_md
intra_log_negMinWP_md_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]/ 
  (meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"] + meanMinWPP50_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"])
mean(intra_log_negMinWP_md_gr.multiMod)
HPDinterval(intra_log_negMinWP_md_gr.multiMod, 0.95)

# log_negP50
intra_log_negP50_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"] + meanMinWPP50_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"])
mean(intra_log_negP50_gr.multiMod)
HPDinterval(intra_log_negP50_gr.multiMod, 0.95)

# Evolutive correlation
cor_log_negMinWP_md_log_negP50_gr.multiMod <- meanMinWPP50_gr.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negP50.animal"]/
  sqrt(meanMinWPP50_gr.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*meanMinWPP50_gr.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(cor_log_negMinWP_md_log_negP50_gr.multiMod)
HPDinterval(cor_log_negMinWP_md_log_negP50_gr.multiMod, 0.95)


############################# HT + HT ~ ENVIRONMENT + group ##############################################################################

# Complete data
plmm.data <- completeFun(plmm.data, c("PC1", "Species"))

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50_env_grNI.multiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:group + (trait:PC1 + trait:PC2 + trait:PC3), 
                                                random = ~us(trait):animal + us(trait):genus.Rand, 
                                                rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                                prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50_env_grNI.multiMod)
plot(log_Ks_log_negP50_env_grNI.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50_env_grNI.multiMod$Sol)
autocorr.diag(log_Ks_log_negP50_env_grNI.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50_env_grNI.multiMod$Sol)
heidel.diag(log_Ks_log_negP50_env_grNI.multiMod$VCV)


## Phylogenetic signal

# log_Ks
physig_log_Ks_env_grNI.multiMod <- log_Ks_log_negP50_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_negP50_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_negP50_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_negP50_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks_env_grNI.multiMod)
HPDinterval(physig_log_Ks_env_grNI.multiMod, 0.95)

# log_negP50
physig_log_negP50_env_grNI.multiMod <- log_Ks_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (log_Ks_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_Ks_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_Ks_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50_env_grNI.multiMod)
HPDinterval(physig_log_negP50_env_grNI.multiMod, 0.95)


## inter-genus

# log_Ks
inter_log_Ks_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(inter_log_Ks_env_gr.multiMod)
HPDinterval(inter_log_Ks_env_gr.multiMod, 0.95)

# log_negP50
inter_log_negP50_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(inter_log_negP50_env_gr.multiMod)
HPDinterval(inter_log_negP50_env_gr.multiMod, 0.95)


## Intra-genus

# log_Ks
intra_log_Ks_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(intra_log_Ks_env_gr.multiMod)
HPDinterval(intra_log_Ks_env_gr.multiMod, 0.95)

# log_negP50
intra_log_negP50_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(intra_log_negP50_env_gr.multiMod)
HPDinterval(intra_log_negP50_env_gr.multiMod, 0.95)

## Evolutive correlation
corr_log_Ks_log_negP50_env_grNI.multiMod <- log_Ks_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_negP50_env_grNI.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_negP50_env_grNI.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(corr_log_Ks_log_negP50_env_grNI.multiMod)
HPDinterval(corr_log_Ks_log_negP50_env_grNI.multiMod , 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_Hv_env_grNI.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:group + (trait:PC1 + trait:PC2 + trait:PC3), 
                                            random = ~us(trait):animal + us(trait):genus.Rand, 
                                            rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                            prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv_env_grNI.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv_env_grNI.multiMod$Sol)
autocorr.diag(log_Ks_log_Hv_env_grNI.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv_env_grNI.multiMod$Sol)
heidel.diag(log_Ks_log_Hv_env_grNI.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks_env_grNI.multiMod <- log_Ks_log_Hv_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_Hv_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_Hv_env_grNI.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks_env_grNI.multiMod)
HPDinterval(physig_log_Ks_env_grNI.multiMod, 0.95)

# log_Hv
physig_log_Hv_env_grNI.multiMod <- log_Ks_log_Hv_env_grNI.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_Ks_log_Hv_env_grNI.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv_env_grNI.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_Ks_log_Hv_env_grNI.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv_env_grNI.multiMod)
HPDinterval(physig_log_Hv_env_grNI.multiMod, 0.95)

# inter-genus

# log_Ks
inter_log_Ks_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(inter_log_Ks_env_gr.multiMod)
HPDinterval(inter_log_Ks_env_gr.multiMod, 0.95)

# log_Hv
inter_log_Hv_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(inter_log_Hv_env_gr.multiMod)
HPDinterval(inter_log_Hv_env_gr.multiMod, 0.95)


# Intra-genus

# log_Ks
intra_log_Ks_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(intra_log_Ks_env_gr.multiMod)
HPDinterval(intra_log_Ks_env_gr.multiMod, 0.95)

# log_Hv
intra_log_Hv_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(intra_log_Hv_env_gr.multiMod)
HPDinterval(intra_log_Hv_env_gr.multiMod, 0.95)

# Evolutive correlation
corr_log_Ks_log_Hv_env_grNI.multiMod  <- log_Ks_log_Hv_env_grNI.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv_env_grNI.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv_env_grNI.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(corr_log_Ks_log_Hv_env_grNI.multiMod)
HPDinterval(corr_log_Ks_log_Hv_env_grNI.multiMod, 0.95)


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_negP50_env_grNI.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1 + trait:group + (trait:PC1 + trait:PC2 + trait:PC3), 
                                                         random = ~us(trait):animal + us(trait):genus.Rand, 
                                                         rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                                         prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_negP50_env_grNI.multiMod)
plot(log_negMinWP_md_log_negP50_env_grNI.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_negP50_env_grNI.multiMod$Sol)
autocorr.diag(log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_negP50_env_grNI.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md_env_grNI.multiMod <- log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md)
HPDinterval(physig_log_negMinWP_md, 0.95)

# log_negP50
physig_log_negP50_env_grNI.multiMod <- log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)

# Evolutive correlation
corr_log_negP50_log_negMinWP_md_env_grNI.multiMod <- log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negP50.animal"]/
  sqrt(log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*log_negMinWP_md_log_negP50_env_grNI.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(corr_log_negP50_log_negMinWP_md_env_grNI.multiMod)
HPDinterval(corr_log_negP50_log_negMinWP_md_env_grNI.multiMod, 0.95)



############################# HT + HT ~ ENVIRONMENT * group ##############################################################################

# Complete data
plmm.data <- completeFun(plmm.data, c("PC1", "Species"))

### log_Ks + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_negP50_env_gr.multiMod <- MCMCglmm(cbind(log_Ks, log_negP50)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                              random = ~us(trait):animal + us(trait):genus.Rand, 
                                              rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                              prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_negP50_env_gr.multiMod)
plot(log_Ks_log_negP50_env_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_negP50_env_gr.multiMod$Sol)
autocorr.diag(log_Ks_log_negP50_env_gr.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_negP50_env_gr.multiMod$Sol)
heidel.diag(log_Ks_log_negP50_env_gr.multiMod$VCV)


## Phylogenetic signal

# log_Ks
physig_log_Ks_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_negP50
physig_log_negP50_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50)
HPDinterval(physig_log_negP50, 0.95)


## Evolutive correlation
cor_log_Ks_log_negP50_env_gr.multiMod <- log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_negP50_env_gr.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_negP50_env_gr.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(corr_log_Ks_log_negP50_env_gr.multiMod)
HPDinterval(corr_log_Ks_log_negP50_env_gr.multiMod, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_Hv_env_gr.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                          random = ~us(trait):animal + us(trait):genus.Rand, 
                                          rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                          prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv_env_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv_env_gr.multiMod$Sol)
autocorr.diag(log_Ks_log_Hv_env_gr.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv_env_gr.multiMod$Sol)
heidel.diag(log_Ks_log_Hv_env_gr.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_Hv
physig_log_Hv_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)


# Evolutive correlation
cor_log_Ks_log_Hv_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv_env_gr.multiMod)
HPDinterval(cor_log_Ks_log_Hv_env_gr.multiMod, 0.95)


### log_negMinWP_md + log_negP50 ###

# Complete data
data <- completeFun(plmm.data, c("log_negP50", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_negP50_env_gr.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_negP50)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                                       random = ~us(trait):animal + us(trait):genus.Rand, 
                                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                                       prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_negP50_env_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_negP50_env_gr.multiMod$Sol)
autocorr.diag(log_negMinWP_md_log_negP50_env_gr.multiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_negP50_env_gr.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_negP50_env_gr.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md_env_gr.multiMod)
HPDinterval(physig_log_negMinWP_md_env_gr.multiMod, 0.95)

# log_negP50
physig_log_negP50_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"]/ 
  (log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(physig_log_negP50_env_gr.multiMod)
HPDinterval(physig_log_negP50_env_gr.multiMod, 0.95)

# inter-genus

# log_negMinWP_md
inter_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"]/ 
  (log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(inter_log_negMinWP_md_env_gr.multiMod)
HPDinterval(inter_log_negMinWP_md_env_gr.multiMod, 0.95)

# log_negP50
inter_log_negP50_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"]/ 
  (log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(inter_log_negP50_env_gr.multiMod)
HPDinterval(inter_log_negP50_env_gr.multiMod, 0.95)

# intra-genus

# log_negMinWP_md
intra_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]/ 
  (log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(intra_log_negMinWP_md_env_gr.multiMod)
HPDinterval(intra_log_negMinWP_md_env_gr.multiMod, 0.95)

# log_negP50
intra_log_negP50_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"]/ 
  (log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.animal"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.genus.Rand"] + log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negP50:traitlog_negP50.units"])
mean(intra_log_negP50_env_gr.multiMod)
HPDinterval(intra_log_negP50_env_gr.multiMod, 0.95)

# Evolutive correlation
cor_log_negP50_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negP50.animal"]/
  sqrt(log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*log_negMinWP_md_log_negP50_env_gr.multiMod$VCV[,"traitlog_negP50:traitlog_negP50.animal"]) 
mean(cor_log_negP50_log_negMinWP_md_env_gr.multiMod)
HPDinterval(cor_log_negP50_log_negMinWP_md_env_gr.multiMod, 0.95)


### log_Ks + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_Ks"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_Ks_log_Hv_env_gr.multiMod <- MCMCglmm(cbind(log_Ks, log_Hv)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                          random = ~us(trait):animal + us(trait):genus.Rand, 
                                          rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                          prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_Ks_log_Hv_env_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_Ks_log_Hv_env_gr.multiMod$Sol)
autocorr.diag(log_Ks_log_Hv_env_gr.multiMod$VCV)
# Convergence
heidel.diag(log_Ks_log_Hv_env_gr.multiMod$Sol)
heidel.diag(log_Ks_log_Hv_env_gr.multiMod$VCV)

# Phylogenetic signal

# log_Ks
physig_log_Ks_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Ks:traitlog_Ks.units"])
mean(physig_log_Ks)
HPDinterval(physig_log_Ks, 0.95)

# log_Hv
physig_log_Hv_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv)
HPDinterval(physig_log_Hv, 0.95)


# Evolutive correlation
cor_log_Ks_log_Hv_env_gr.multiMod <- log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Ks.animal"]/
  sqrt(log_Ks_log_Hv_env_gr.multiMod$VCV[, "traitlog_Ks:traitlog_Ks.animal"]*log_Ks_log_Hv_env_gr.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Ks_log_Hv_env_gr.multiMod)
HPDinterval(cor_log_Ks_log_Hv_env_gr.multiMod, 0.95)


### log_negMinWP_md + log_Hv ###

# Complete data
data <- completeFun(plmm.data, c("log_Hv", "log_negMinWP_md"))
# Intersenction data
phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% data$genus])
data <- data[data$genus %in% phylo$tip.label, ]
length(data$genus)

# Model
log_negMinWP_md_log_Hv_env_gr.multiMod <- MCMCglmm(cbind(log_negMinWP_md, log_Hv)~trait-1 + trait:group * (trait:PC1 + trait:PC2 + trait:PC3), 
                                                       random = ~us(trait):animal + us(trait):genus.Rand, 
                                                       rcov = ~us(trait):units, pedigree=phylo, family = c("gaussian","gaussian"),
                                                       prior = priorMulti, data= data, nitt = nitter, burnin = burni, thin = thinni)

summary(log_negMinWP_md_log_Hv_env_gr.multiMod)

## Model diagnostics

# Autocorrelation
autocorr.diag(log_negMinWP_md_log_Hv_env_gr.multiMod$Sol)
autocorr.diag(log_negMinWP_md_log_Hv_env_gr.multiMod$VCV)
# Convergence
heidel.diag(log_negMinWP_md_log_Hv_env_gr.multiMod$Sol)
heidel.diag(log_negMinWP_md_log_Hv_env_gr.multiMod$VCV)


# Phylogenetic signal

# log_negMinWP_md
physig_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]/ 
  (log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(physig_log_negMinWP_md_env_gr.multiMod)
HPDinterval(physig_log_negMinWP_md_env_gr.multiMod, 0.95)

# log_Hv
physig_log_Hv_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"]/ 
  (log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(physig_log_Hv_env_gr.multiMod)
HPDinterval(physig_log_Hv_env_gr.multiMod, 0.95)

# inter-genus

# log_negMinWP_md
inter_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"]/ 
  (log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(inter_log_negMinWP_md_env_gr.multiMod)
HPDinterval(inter_log_negMinWP_md_env_gr.multiMod, 0.95)

# log_Hv
inter_log_Hv_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"]/ 
  (log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(inter_log_Hv_env_gr.multiMod)
HPDinterval(inter_log_Hv_env_gr.multiMod, 0.95)

# intra-genus

# log_negMinWP_md
intra_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"]/ 
  (log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.animal"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.genus.Rand"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_negMinWP_md:traitlog_negMinWP_md.units"])
mean(intra_log_negMinWP_md_env_gr.multiMod)
HPDinterval(intra_log_negMinWP_md_env_gr.multiMod, 0.95)

# log_Hv
intra_log_Hv_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"]/ 
  (log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.animal"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.genus.Rand"] + log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_Hv:traitlog_Hv.units"])
mean(intra_log_Hv_env_gr.multiMod)
HPDinterval(intra_log_Hv_env_gr.multiMod, 0.95)

# Evolutive correlation
cor_log_Hv_log_negMinWP_md_env_gr.multiMod <- log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_Hv.animal"]/
  sqrt(log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[, "traitlog_negMinWP_md:traitlog_negMinWP_md.animal"]*log_negMinWP_md_log_Hv_env_gr.multiMod$VCV[,"traitlog_Hv:traitlog_Hv.animal"]) 
mean(cor_log_Hv_log_negMinWP_md_env_gr.multiMod)
HPDinterval(cor_log_Hv_log_negMinWP_md_env_gr.multiMod, 0.95)


#### SAVE IMAGE ####
save.image(file = "RData/Phylogenetic_linear_mixed_multi_models_groups.RData")
