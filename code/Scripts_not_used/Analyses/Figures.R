############################################# MCMC REGRESSION PLOTS ########################################################################

#### PACKAGES ####
install.packages("MCMCglmmExtras", dependencies = T)
library(ggbiplot)
library(ggplot2)
library(factoextra)
library(MCMCglmmExtras)

library(dplyr)
library(plyr)
library(bayr)
library(postMCMCglmm)
library(ggsci)
library(scales)

library(ggmcmc)
library(cowplot)

setwd("~/Pablo/data")

# Palette 
show_col(pal_npg("nrc")(10))

#### R DATA ####

## PCA
load("Data/RData/PCA.RData")

## MCMCGLMM
#Non phylogenetic
load("RData/Non_Phylogenetic_linear_mixed_multi_models_groups.RData")
# Whitout group
load("RData/Phylogenetic_linear_mixed_multi_models.RData")
# With group
load("RData/Phylogenetic_linear_mixed_multi_models_groups.RData")


#### FUNCTIONS ####

## PCAs
pcaPlot <- function(pca, axs){
p <- fviz_pca_biplot(pca, axes = axs, label = "none",
                     geom.ind = "point", 
                     habillage = ht_pca_gr, alpha.ind = 0.5,
                     addEllipses=TRUE, ellipse.level=0.95,
                     # select.var = list(contrib = 10),
                     col.var = "#000000", 
                     repel = T,
                     title = "",
                     theme(axis.title.x =element_blank(), axis.title.y= element_blank())
)
p <- p + xlim(-8, 8) + ylim (-8, 8)
p <- p + scale_color_npg() +  theme_minimal() + theme(legend.position = "none")
print(p)
}

# PCA-Ht correlations
PcaCorPlot <- function(pca, axs){
  p <- fviz_pca_ind(pca, geom = "point", col.ind = ht_pca_de, alpha.ind = 0.5,
                    addEllipses=TRUE, ellipse.level=0.95, title = ""
  )
p <- p + scale_color_npg() +  theme_minimal() + theme(legend.position = "none")
p <- fviz_add(p, hyd.coord[, axs]*10, addlabel = F, 
              geom="arrow", repel = T, linetype = "solid", color = "darkblue", title = "")
p <- p + xlim(-8, 8) + ylim (-8, 8)
print(p)
}

## mcmcReg Data preparation
mcmcRegDataFun <- function(mdl, vars){
  # Dataframes
  postPlot <- data.frame( mean = numeric(), low95CI = numeric(), high95CI = numeric(), sign = factor(), stringsAsFactors=FALSE)
  # Intercept
  mean <- mean(mdl$model$Sol[,"(Intercept)"])
  CI <- HPDinterval(mdl$model$Sol[,"(Intercept)"], prob=0.95)
  if(CI[1] < 0 && CI[2] > 0){
    sign <- "No"
  } else {
    sign <- "Yes"
  }
  post.df <- data.frame(
    mean = mean, 
    low95CI = CI[1], 
    high95CI = CI[2],
    sign = sign)
  row.names(post.df) <- "Intercept"
  postPlot <- rbind(postPlot, post.df)
  # Fixed effects
  for(var in vars){
    mean <- mean(mdl$model$Sol[, var])
    CI <- HPDinterval(mdl$model$Sol[,var],prob=0.95)
    if(CI[1] < 0 && CI[2] > 0){
      sign <- "No"
    } else {
      sign <- "Yes"
    }
    post.df <- data.frame(
      mean = mean, 
      low95CI = CI[1], 
      high95CI = CI[2],
      sign = sign)
    row.names(post.df) <- var
    postPlot <- rbind(postPlot, post.df)
  } # end of fixed effects
  return(postPlot)
} # end of function

## Univariate plot
mcmcRegPlotFun <- function(pointsData = plmm.data, mdl, mdl_gr, NoPhy_mdl, varX, varY, CI = F, xmin = -8, xmax = 8, ylim, ymin = -4, ymax = 4){
  # Model selection
  dic <- DIC(mdl$model, mdl_gr$model)
  deltaDic <- dic[1, "DIC"] - dic[2, "DIC"] 
  # No groups model
  if(!isTRUE((deltaDic) >= 4)){
    # Coefficients
    coefData <- mcmcRegDataFun(mdl, vars = c(varX))
    data <- completeFun(pointsData, c(varY, varX))
    # Plot
    p <- ggplot(aes(x = data[, varX], y = data[, varY], col = group), data = data) + geom_point(size = 1.2, alpha = 0.4)
    p <- p + geom_abline(aes(intercept = coefData["Intercept", "mean"], slope = coefData[varX, "mean"]), 
                         lty = 1, lwd = 1.2, data = coefData, color = "black")
    # p <- p + xlab(paste(varX)) + ylab(paste(varY))
    # Mean CI
    if(CI == TRUE){
    p <- p + geom_abline(aes(intercept = coefData["Intercept", "low95CI"], slope = coefData[varX, "mean"]),
                             lty = 4, lwd = 0.8, data = coefData, color = "black")
    p <- p + geom_abline(aes(intercept = coefData["Intercept", "high95CI"], slope = coefData[varX, "mean"]),
                         lty = 4, lwd = 0.8, data = coefData, color = "black")
    }
  } else {
    # Group models
    # Coefficients
    varXGym <- paste0(varX, ":groupGymnosperms")
    coefData_gr <- mcmcRegDataFun(mdl_gr, vars = c(varX, "groupGymnosperms", varXGym))
    gimInt <- coefData_gr["Intercept", "mean"] + coefData_gr["groupGymnosperms", "mean"]
    gymSlope <- coefData_gr[varX, "mean"] + coefData_gr[varXGym, "mean"]
    data <- completeFun(pointsData, c(varY, varX))
    # Angiosperms
    p <- ggplot(aes(x = data[, varX], y = data[, varY], col = group), data = data) + geom_point(size = 1.2, alpha = 0.4)
    p <- p + geom_abline(aes(intercept = coefData_gr["Intercept", "mean"], slope = coefData_gr[varX, "mean"]), 
                         lty = 1, lwd = 1.2, data = coefData_gr, color = "#00A087FF")
    p <- p + xlab(paste(varX)) + ylab(paste(varY))
    # # Mean CI
    if(CI == TRUE){
      p <- p + geom_abline(aes(intercept = coefData_gr["Intercept", "low95CI"], slope = coefData_gr[varX, "mean"]),
                           lty = 4, lwd = 0.8, data = coefData_gr, color = "#00A087FF")
      p <- p + geom_abline(aes(intercept = coefData_gr["Intercept", "high95CI"], slope = coefData_gr[varX, "mean"]),
                           lty = 4, lwd = 0.8, data = coefData_gr, color = "#00A087FF")
    }
    # Gymnosperms slope
    if(coefData_gr[varXGym, "sign"] == "Yes"){
      p <- p + geom_abline(aes(intercept = gimInt, 
                               slope = gymSlope), col = "#3C5488FF", lty = 1, lwd = 1, data = coefData_gr)
      # Mean CI
      if(CI == TRUE){
        gymLowCI <- coefData_gr["Intercept", "low95CI"] + coefData_gr["groupGymnosperms", "low95CI"]
        gymHighCI <- coefData_gr["Intercept", "high95CI"] + coefData_gr["groupGymnosperms", "high95CI"]
        p <- p + geom_abline(aes(intercept = gymLowCI, slope = gymSlope),
                             lty = 4, lwd = 0.8, data = coefData_gr, color = "#3C5488FF")
        p <- p + geom_abline(aes(intercept = gymHighCI, slope = gymSlope),
                             lty = 4, lwd = 0.8, data = coefData_gr, color = "#3C5488FF")
      }
    } else {
    # Gymnosperm intercept
    if(coefData_gr["groupGymnosperms", "sign"] == "Yes"){
      p <- p + geom_abline(aes(intercept = gimInt,
                               slope = coefData_gr[varX, "mean"]), data = coefData_gr,
                           color="#3C5488FF", size = 1.2)
    } # End only gymnosperm intercept
   } # End model with groups
  }
  p <- p + scale_color_npg() + theme_minimal() + 
    theme(legend.position = "none") + theme( axis.line = element_line(colour = "darkgrey", 
                                                                      size = 1, linetype = "solid"))
  p <- p + theme(axis.title.x =element_blank(), axis.title.y= element_blank()) + xlim(xmin, xmax)
  if(isTRUE(ylim)){ 
  p <- p + ylim(ymin, ymax)
  }
  # Non phylogenetic regression
  noPhyData <-  mcmcRegDataFun(NoPhy_mdl, vars = c(varX))
  p <- p + geom_abline(aes(intercept = noPhyData["Intercept", "mean"], slope = noPhyData[varX, "mean"]),
                       lty = "dashed", lwd = 0.8, data = noPhyData, color = "#000000")
  plot(p, res = 1000)
}

## Correlation densities
corDensFun <- function(chain, chain2, name_cor){
  # Dataframe
  df1 <- data.frame("Iteration" = 1:length(chain), 
                    "Chain" = as.character(rep(1, length = length(chain))),
                    "Parameter" = rep("cor", length = length(chain)),
                    "var1" = chain)
  df2 <- data.frame("Iteration" = 1:length(chain2), 
                    "Chain" = as.character(rep(2, length = length(chain2))),
                    "Parameter" = rep("cor_env", length = length(chain2)),
                    "value" = chain2)
  df <- bind_rows(df1, df2)
  names(df)[names(df) == "var1"] <- "Cor"
  # Mean and HDP intervals
  mu <- ddply(df, "Parameter", summarise, mean = mean(Cor))
  mu[1, "low_Hdp"] <- HPDinterval(chain)[1]
  mu[1, "high_Hdp"] <- HPDinterval(chain)[2]
  mu[2, "low_Hdp"] <- HPDinterval(chain2)[1]
  mu[2, "high_Hdp"] <- HPDinterval(chain2)[2]
  # Density plot
  p<-ggplot(data = df, aes(x= Cor, fill = Parameter, col = Parameter)) + geom_density(alpha = 0.4)+
    geom_vline(data = mu, aes(xintercept = mean, color = Parameter), linetype="dashed", lwd = 0.8) + 
    geom_vline(data = mu, aes(xintercept = low_Hdp, color = Parameter), linetype= 3, lwd = 0.8) + 
    geom_vline(data = mu, aes(xintercept = high_Hdp, color = Parameter), linetype= 3, lwd = 0.8) +
    geom_vline(aes(xintercept = 0), col = "black", linetype= 1, lwd = 0.5)
  p <- p + scale_color_npg() + theme_minimal() + theme(legend.position="none") + xlab(name_cor)
  plot(p, res = 1000)
}

## Trace plot 
corTraceFun <- function(chain, chain2, name_cor){
  df1 <- data.frame("Iteration" = 1:length(chain), 
                    "Chain" = as.character(rep(1, length = length(chain))),
                    "Parameter" = rep("cor", length = length(chain)),
                    "var1" = chain)
  df2 <- data.frame("Iteration" = 1:length(chain2), 
                    "Chain" = as.character(rep(2, length = length(chain2))),
                    "Parameter" = rep("cor_env", length = length(chain2)),
                    "value" = chain2)
  df <- bind_rows(df1, df2)
  names(df)[names(df) == "var1"] <- "Cor"
  # Trace Plot
  p <- ggplot(df, aes(x=Iteration, y=Cor, colour=as.factor(Chain)))
  p <- p + geom_line(alpha=0.5) + scale_colour_discrete(name="Chain")
  p <- p + scale_color_npg() +  theme_minimal() + theme(legend.position="none") + xlab(name_cor) +
  geom_hline(aes(yintercept = 0), col = "black", linetype= 1, lwd = 0.5)
  plot(p, res = 1000)
}

## Effects means
plot.estimates <- function(mdl, gr, var1, var2, fixed, random, resid, interaction) {
  x <- summary(mdl)
  if(isTRUE(gr)){
  # Fixed (interaction)
  df.fix <- rbind(x$solutions[1:10, 1:4])
  rownames(df.fix) <- c(paste0(var1, "Ang"), paste0(var2, "Ang"), paste0(var1, "Gym"), paste0(var2, "Gym"), paste0(var1, ":PC1"), paste0(var2, ":PC1"), 
                        paste0(var1, ":PC2"), paste0(var2, ":PC2"), paste0(var1, ":PC3"), paste0(var2, ":PC3"))
  } else {
    df.fix <- rbind(x$solutions[1:8, 1:4])
    rownames(df.fix) <- c(paste0(var1), paste0(var2), paste0(var1, ":PC1"), paste0(var2, ":PC1"), 
                          paste0(var1, ":PC2"), paste0(var2, ":PC2"), paste0(var1, ":PC3"), paste0(var2, ":PC3"))
  }
  # Interaction
  if(interaction == T){
    df.fixInter <- rbind(x$solutions[, 1:4])
    rownames(df.fixInter) <- c(paste0(var1, "Ang"), paste0(var2, "Ang"), paste0(var1, ":Gym"), paste0(var2, ":Gym"), 
                               paste0(var1, ":PC1"), paste0(var2, ":PC1"), paste0(var1, ":PC2"), paste0(var2, ":PC2"), paste0(var1, ":PC3"),
                               paste0(var2, ":PC3"), paste0(var1, ":Gym:PC1"), paste0(var2, ":Gym:PC1"), paste0(var1, ":Gym:PC2"), 
                               paste0(var2, ":Gym:PC2"), paste0(var1, ":Gym:PC3"), paste0(var2, ":Gym:PC3"))
  }
  # Phylogeny
  df.rand <- rbind(x$Gcovariances[c(1, 2, 4),])
  rownames(df.rand) <- c(paste0(var1, "Var"), paste0(var2, ":", var1, "Cov"), paste0(var2, "Var"))
  df.rand <- df.rand[c(paste0(var1, "Var"), paste0(var2, "Var"), paste0(var2, ":", var1, "Cov")),]
  # Residuals 
  df.res <- rbind(x$Rcovariances)
  rownames(df.res) <- c(paste0(var1, "Res"), paste0(var2, ":", var1, "Res"), paste0(var1, ":", var2, "Res"), paste0(var2, "Res"))
  # Fixed effects
  if(fixed == TRUE){
    df <- df.fix
    # Fixed interactions
    if(interaction == TRUE){
      df <- df.fixInter
    }
  }  # end fix
  # Random effects only
  if(random == TRUE){
    df <- df.rand
  }
  # fixed + random
  if(fixed == TRUE && random == TRUE)
    if(interaction == TRUE) {
      df <- rbind(df.fixInter, df.rand)
    } else {
      df <- rbind(df.fix, df.rand)
    }
  # Fixed + random + residuals
  if(fixed == TRUE && random == TRUE && resid == TRUE){
    if(interaction == TRUE) {
      df <- rbind(df.fixInter, df.rand, df.res)
    } else {
      df <- rbind(df.fix, df.rand, df.res)
    }
  }
  df <- as.data.frame(df)
  # Significanace
  for(i in 1:length(df[, 1])){
    if(isTRUE(df[i, 2] <= 0 && df[i, 3] > 0)){ 
      df[i, "sgnf"] <- "No"
  } else {
      df[i, "sgnf"] <- "Yes"
  }
}
  sgnf <- as.factor(df$sgnf)
  n <- dim(df)[1]
  par(mar=c(2, 10, 4, 1))
  plot(df[,1], n:1, yaxt="n", ylab="", xlim=range(-2.5, 2.5), pch=19,   #xlim=range(min(df[, 2]) * 1.2, max(df[, 3]) * 1.2)
       col = c("black", "#DC0000FF")[sgnf], grid())
  axis(2, at=n:1, labels = rownames(df), las=2)
  arrows(df[,2], n:1, df[,3], n:1, code=0, col = c("black", "#DC0000FF")[sgnf])
  abline(v=0, lty=2)
  if(random == T){
    abline(h = 4.5, lty = 4, col = "#DC0000FF")
  }
  if(random == T && resid == T){
    abline(h = 8.5, lty = 4, col = "#DC0000FF")  
  }
}

#### PCA ####

## Biplot

# PC1 -  PC2
f.out <- paste0("C:/Users/pablo/Desktop/Master/TFM/Analisis/results/pca/Env_pca_Hydraulic_groups.pdf")
pdf(f.out, f.out, width=5, height=5)
pcaPlot(pca, axs = c(1,2))
dev.off()
print(paste("==>", f.out))

# PC1- PC3
f.out <- paste0("C:/Users/pablo/Desktop/Master/TFM/Analisis/results/pca/Env_pca_Hydraulic_groups_pc3.pdf")
pdf(f.out, f.out, width=5, height=5)
pcaPlot(pca, axs = c(1,3))
dev.off()
print(paste("==>", f.out))

# Hydraulic traits PC1-2
f.out <- paste0("C:/Users/pablo/Desktop/Master/TFM/Analisis/results/pca/Hyd_corr_pc1_pc2.pdf")
pdf(f.out, f.out, width=5, height=5)
PcaCorPlot(pca, axs = c(1, 2))
dev.off()
print(paste("==>", f.out))

# Hydraulic traits PC1-3
f.out <- paste0("C:/Users/pablo/Desktop/Master/TFM/Analisis/results/pca/Hyd_corr_pc1_pc3.pdf")
pdf(f.out, f.out, width=5, height=5)
PcaCorPlot(pca, axs = c(1, 3))
dev.off()
print(paste("==>", f.out))


### LEGEND ####

f.out <- paste0("C:/Users/pablo/Desktop/Master/TFM/Analisis/results/AngGym_legend.pdf")
pdf(f.out, f.out, width=10, height=10)
legend <- get_legend(p)
plot(legend, res = 1000)
dev.off()
print(paste("==>", f.out))


#### Ht ~ Env Plots ####

### Log_Ks

f.out <- paste0("results/MCMCreg/log_Ks_env/log_Ks_PC1_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Ks_PC1.mod, mdl_gr = log_Ks_PC1_gr.mod,
               varX = "PC1", varY = "log_Ks",  CI = F, ylim = F)
paste("======> Image saved on ", f.out)
dev.off()


f.out <- paste0("results/MCMCreg/log_Ks_env/log_Ks_PC2_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Ks_PC2.mod, mdl_gr = log_Ks_PC2_gr.mod, 
               varX = "PC2", varY = "log_Ks", CI = F, ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_Ks_env/log_Ks_PC3_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Ks_PC3.mod, mdl_gr = log_Ks_PC3_gr.mod, 
               varX = "PC3", varY = "log_Ks", CI = F, ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

### Log_negP50

f.out <- paste0("results/MCMCreg/log_negP50_env/log_negP50_PC1_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negP50_PC1.mod, mdl_gr = log_negP50_PC1_gr.mod, 
               varX = "PC1", varY = "log_negP50", CI = F, ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_negP50_env/log_negP50_PC2_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negP50_PC2.mod, mdl_gr = log_negP50_PC2_gr.mod, 
               varX = "PC2", varY = "log_negP50", CI = F, ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_negP50_env/log_negP50_PC3_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negP50_PC3.mod, mdl_gr = log_negP50_PC3_gr.mod, 
               varX = "PC3", varY = "log_negP50", CI = F, ylim = F)
paste("======> Image saved on ", f.out)
dev.off()


### log_negMinWP_md

f.out <- paste0("results/MCMCreg/log_negMinWP_md_env/log_negMinWP_md_PC1_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negMinWP_md_PC1.mod, mdl_gr = log_negMinWP_md_PC1_gr.mod,
               varX = "PC1", varY = "log_negMinWP_md", ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_negMinWP_md_env/log_negMinWP_md_PC2_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negMinWP_md_PC2.mod, mdl_gr = log_negMinWP_md_PC2_gr.mod,
               varX = "PC2", varY = "log_negMinWP_md", ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_negMinWP_md_env/log_negMinWP_md_PC3_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negMinWP_md_PC3.mod, mdl_gr = log_negMinWP_md_PC3_gr.mod, 
               varX = "PC3", varY = "log_negMinWP_md", ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

### log_Hv

f.out <- paste0("results/MCMCreg/log_Hv_env/log_Hv_PC1_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Hv_PC1.mod, mdl_gr = log_Hv_PC1_gr.mod, 
               varX = "PC1", varY = "log_Hv", ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_Hv_env/log_Hv_PC2_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Hv_PC2.mod, mdl_gr = log_Hv_PC2_gr.mod, 
               varX = "PC2", varY = "log_Hv", ylim = F)
paste("======> Image saved on ", f.out)
dev.off()

f.out <- paste0("results/MCMCreg/log_Hv_env/log_Hv_PC3_MCMCreg.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Hv_PC3.mod, mdl_gr = log_Hv_PC3_gr.mod, 
               varX = "PC3", varY = "log_Hv", ylim = F)
paste("======> Image saved on ", f.out)
dev.off()


#### Ht ~ Ht ####

## log_Ks ~ log_negP50

f.out <- paste0("results/MCMCreg/trait_trait/log_Ks_log_negP50.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Ks_log_negP50.mod, mdl_gr = log_Ks_log_negP50_gr.mod, NoPhy_mdl = log_Ks_log_negP50_NoPhy.mod,
               varY = "log_Ks", varX = "log_negP50", xmin = -5, xmax = 5, ylim = T,  ymin = -5, ymax = 5)
paste("======> Image saved on ", f.out)
dev.off()

## log_negP50 ~ log_Ks

f.out <- paste0("results/MCMCreg/trait_trait/log_negP50_log_Ks.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negP50_log_Ks.mod, mdl_gr = log_negP50_log_Ks_gr.mod, NoPhy_mdl = log_negP50_log_Ks_NoPhy.mod,
               varY = "log_negP50", varX = "log_Ks", xmin = -5, xmax = 5, ylim = T,  ymin = -5, ymax = 5)
paste("======> Image saved on ", f.out)
dev.off()


## log_Hv ~ log_Ks

f.out <- paste0("results/MCMCreg/trait_trait/log_Hv_log_Ks.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Hv_log_Ks.mod, mdl_gr = log_Hv_log_Ks_gr.mod, NoPhy_mdl = log_Hv_log_Ks_NoPhy.mod,
               varY = "log_Hv", varX = "log_Ks", xmin = -5, xmax = 5, ylim = T,  ymin = -5, ymax = 5)
paste("======> Image saved on ", f.out)
dev.off()

## log_Ks ~ log_Hv

f.out <- paste0("results/MCMCreg/trait_trait/log_Ks_log_Hv.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_Ks_log_Hv.mod, mdl_gr = log_Ks_log_Hv_gr.mod, NoPhy_mdl = log_Ks_log_Hv_NoPhy.mod,
               varY = "log_Ks", varX = "log_Hv", xmin = -5, xmax = 5, ylim = T,  ymin = -5, ymax = 5)
paste("======> Image saved on ", f.out)
dev.off()


## log_negP50 ~ log_negMinWP_md

f.out <- paste0("results/MCMCreg/trait_trait/log_negP50_log_negMinWP_md.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negP50_log_negMinWP_md.mod, mdl_gr = log_negP50_log_negMinWP_md_gr.mod, NoPhy_mdl = log_negP50_log_negMinWP_md_NoPhy.mod,
               varY = "log_negP50", varX = "log_negMinWP_md", xmin = -5, xmax = 5, ylim = T,  ymin = -3, ymax = 3)
paste("======> Image saved on ", f.out)
dev.off()

## log_negMinWP_md ~ log_negP50

f.out <- paste0("results/MCMCreg/trait_trait/log_negMinWP_md_log_negP50.pdf")
pdf(f.out, width=2.5, height=2.5)
mcmcRegPlotFun(pointsData = plmm.data, mdl = log_negMinWP_md_log_negP50.mod, mdl_gr = log_negMinWP_md_log_negP50_gr.mod, NoPhy_mdl = log_negMinWP_md_log_negP50_NoPhy.mod,
               varY = "log_negMinWP_md", varX = "log_negP50", xmin = -5, xmax = 5, ylim = T,  ymin = -3, ymax = 3)
paste("======> Image saved on ", f.out)
dev.off()

### Correlation plots ####

### Ht - Ht densities

## log_Ks - log_negP50
f.out <- paste0("results/Evo_corr/dens/log_Ks_log_negP50_Dens.pdf")
pdf(f.out, width=3, height=3)
corDensFun(chain = cor_log_Ks_log_negP50.multiMod, chain2 = cor_log_Ks_log_negP50_env.multiMod, 
           name_cor = "Correlation")
paste("======> Image saved on ", f.out)
dev.off()

## log_Ks - log_Hv
f.out <- paste0("results/Evo_corr/dens/log_Ks_log_Hv_Dens.pdf")
pdf(f.out, width=3, height=3)
corDensFun(chain = cor_log_Ks_log_Hv.multiMod, chain2 = cor_log_Ks_log_Hv_env.multiMod, 
           name_cor = "Correlation")
paste("======> Image saved on ", f.out)
dev.off()

## log_minWP_md - log_negP50
f.out <- paste0("results/Evo_corr/dens/log_minWP_md_log_negP50_Dens.pdf")
pdf(f.out, width=3, height=3)
corDensFun(chain = cor_log_negMinWP_md_log_negP50_gr.multiMod, chain2 = cor_log_negP50_log_negMinWP_md_env_gr.multiMod, 
           name_cor = "Correlation")
paste("======> Image saved on ", f.out)
dev.off()


### Correlation trace plots

## log_Ks - log_negP50
f.out <- paste0("results/Evo_corr/trace/log_Ks_log_negP50_Trace.pdf")
pdf(f.out, width=3, height=3)
corTraceFun(chain = cor_log_Ks_log_negP50.multiMod, chain2 = cor_log_Ks_log_negP50_env.multiMod, 
           name_cor = "Iterations")
paste("======> Image saved on ", f.out)
dev.off()

## log_Ks - log_Hv
f.out <- paste0("results/Evo_corr/trace/log_Ks_log_Hv_Trace.pdf")
pdf(f.out, width=3, height=3)
corTraceFun(chain = cor_log_Ks_log_Hv.multiMod, chain2 = cor_log_Ks_log_Hv_env.multiMod, 
            name_cor = "Iterations")
paste("======> Image saved on ", f.out)
dev.off()

## log_minWP_md - log_negP50
f.out <- paste0("results/Evo_corr/trace/log_minWP_log_negP50_Trace.pdf")
pdf(f.out, width=3, height=3)
corTraceFun(chain = cor_log_negMinWP_md_log_negP50_gr.multiMod, chain2 = cor_log_negP50_log_negMinWP_md_env_gr.multiMod, 
            name_cor = "Iterations")
paste("======> Image saved on ", f.out)
dev.off()


### Effects means ####

## log_Ks, log_negP50 fixed effects (no group)

f.out <- paste0("results/MCMCreg/multi_trait_effects/log_Ks_log_negP50_effects.pdf")
pdf(f.out, width=5, height=5)
plot.estimates(log_Ks_log_negP50.multiMod, gr = F, var1 = "log(Ks)", var2= "log(-P50)", fixed = T, random = T, resid = F, interaction = F)
paste("======> Image saved on ", f.out)
dev.off()

## log_Ks, log_Hv fixed effects (no group)

f.out <- paste0("results/MCMCreg/multi_trait_effects/log_Ks_log_Hv_effects.pdf")
pdf(f.out, width=5, height=5)
plot.estimates(log_Ks_log_Hv.multiMod, gr = F, var1 = "log(Ks)", var2= "log(Hv)", fixed = T, random = T, resid = F, interaction = F)
paste("======> Image saved on ", f.out)
dev.off()

## log_negMinWP_md, log_negP50 fixed effects

f.out <- paste0("results/MCMCreg/multi_trait_effects/log_negMinWP_md_log_negP50_effects.pdf")
pdf(f.out, width=5, height=5)
plot.estimates(log_negMinWP_md_log_negP50_env_gr.multiMod, gr = T, var1 = "log(-minWP)", var2= "log(-P50)",
               fixed = T, random = T, resid = F, interaction = T)
paste("======> Image saved on ", f.out)
dev.off()

#### SUPPLEMENTARIES ####

### PPCA ####
load("RData/Phylogenetic_pca.RData")

biplot(ppca, pch = 19)


fviz_pca_biplot(pca, axes = c(1,2), label = "none",
                geom.ind = "point", 
                habillage = ht_pca_gr, alpha.ind = 0.5,
                addEllipses=TRUE, ellipse.level=0.95,
                # select.var = list(contrib = 10),
                col.var = "#000000", 
                repel = T,
                title = "",
                theme(axis.title.x =element_blank(), axis.title.y= element_blank())
)
p <- p + xlim(-8, 8) + ylim (-8, 8)
p <- p + scale_color_npg() +  theme_minimal() + theme(legend.position = "none")

