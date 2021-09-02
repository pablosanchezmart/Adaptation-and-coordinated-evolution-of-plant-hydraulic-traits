##### VARIABLE TRANSFORMATIONS #############################################################################

# Pablo Sanchez-Martinez

print("Checking variables normality and transforming some of them to approach normality")

### PACKAGES ####
# devtools::install_github("kassambara/ggpubr")

library(ggplot2)
library(ggpubr)
library(dplyr)
library(Hmisc)

library(nortest)

library(gam)
library(knitr)

### FUNCTIONS ----------------------------------------------------------------------------------------- ####

transfPlotFun <- function(df = ht, var = "MAT"){
  dens <- ggdensity(df[, var], main = paste0("Density plot of ", var), xlab = var)
  qq <- ggqqplot(df[, var], main = paste0("QQ plot of ", var))
  
  ## Log transformation
  if(all(na.omit(df[, var] < 0))){
    df[, var] <- -df[, var]
    }
  log.vec <- log(df[, var])
  dens_log <- ggdensity(log.vec, main = paste0("Density plot of log(", var, ")"), xlab = paste0("log(", var, ")"))
  qq_log <- ggqqplot(log.vec, main = paste0("QQ plot of log(", var, ")"))
  
  ## Sqrt transformation
  sqrt.vec <- sqrt(df[, var])
  dens_sqrt <- ggdensity(sqrt.vec, main = paste0("Density plot of sqrt(", var, ")"), xlab = paste0("sqrt(", var, ")"))
  qq_sqrt <- ggqqplot(sqrt.vec, main = paste0("QQ plot of sqrt(", var, ")"))
  
  p <- ggarrange(dens, qq, dens_log, qq_log, dens_sqrt, qq_sqrt)
  print(p)
}

### DATA ---------------------------------------------------------------------------------------------- ####

ht <- read.csv("data/HydraTRY_Pablo_2020.csv", header = T)

#### HYDRAULIC TRAITS --------------------------------------------------------------------------------- ####

# P50
transfPlotFun(df = ht, var = "P50")

# MinWP_md
transfPlotFun(df = ht, var = "MinWP_md")

# Ks
transfPlotFun(df = ht, var = "Ks")

# Hv
transfPlotFun(df = ht, var = "Hv")

# Kl
transfPlotFun(df = ht, var = "Kl")

# HSM
transfPlotFun(df = ht, var = "HSM")

#### ENVIRONMENTAL VARS ----------------------------------------------------------------------------------- ####

names(ht)
vars <- names(ht)[35:68]
for(var in vars){
  print(var)
  transfPlotFun(df = ht, var = var)
}

#### VARIABLE TRANSFORMATIONS --------------------------------------------------------------------------- ####

ht$log_negP50 <- log(-ht$P50)
ht$log_negMinWP_md <- log(-ht$MinWP_md)
ht$log_Ks <- log(ht$Ks)
ht$log_Hv <- log(ht$Hv)
ht$log_Kl <- log(ht$Kl)

ht$log_TRange <- log(ht$TRange)
ht$log_TSeasonality <- log(ht$TSeasonality)
ht$sqrt_Prec_ColdQ <- sqrt(ht$Prec_ColdQ)
ht$sqrt_Prec_WarmQ <- sqrt(ht$Prec_WarmQ)
ht$sqrt_Prec_DryQ <- sqrt(ht$Prec_DryQ)
ht$log_Prec_WetQ <- log(ht$Prec_WetQ)
ht$log_Prec_Seasonality <- log(ht$Prec_Seasonality)
ht$sqrt_Prec_Driest <- sqrt(ht$Prec_Driest)
ht$log_Prec_Wettest <- log(ht$Prec_Wettest)

# Factors
names(ht)
facVars <- names(ht)[1:9]
# Factors
for(i in facVars){
  ht[, i] <- as.factor(ht[, i])
}

# Save final dataset

write.csv(ht, "data/HydraTRY_Pablo_2020.csv", row.names = F)
names(ht)
length(ht$Species)
