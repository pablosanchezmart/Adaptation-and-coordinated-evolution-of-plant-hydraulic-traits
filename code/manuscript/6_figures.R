#### MANUSCRIPT FIGURES ###################################################################################################################

# P Sanchez-Martinez

print("Preparing figures...")
remove(list = ls())

# RData

print("Loading all models results...")
load("data/RData/plmm_univariate_NP_results_diagn.RData")
load("data/RData/plmm_univariate_results_diagn.RData")

load("data/RData/plmm_multivariate_results_diagn.RData")
load("data/RData/plmm_multivariate_NP_results_diagn.RData")

load("data/RData/plmm_univariate_Spp_results_diagn.RData")
load("data/RData/plmm_multivariate_Spp_results_diagn.RData")

load("data/RData/PCA.RData")

source("code/manuscript/1_init.R")

### Palette
show_col(pal_npg("nrc")(10))

#### FUNCTIONS --------------------------------------------------------------------------------------------------------####

## Auxiliar functions 

# Complete cases
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
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

## Plot functions

# PCAs
pcaPlot <- function(pca, axs){
p <- fviz_pca_biplot(pca, axes = axs, label = "none",
                     geom.ind = "point", 
                     habillage = ht_pca_gr, alpha.ind = 0.5,
                     addEllipses=TRUE, ellipse.level=0.95,
                     select.var = list(contrib = 10),
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

# mcmcReg Data preparation
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

# Univariate plot
mcmcRegPlotFun <- function(pointsData = plmm.data, xmin = -8, xmax = 8, ylim = T, ymin = -10, ymax = 10){
  # No groups model
  if(!is.na(str_extract(mdl.name, "NG"))){
    # Coefficients
    coefData <- mcmcRegDataFun(mod, vars = c(varX))
    data <- completeFun(pointsData, c(varY, varX))
    # Plot
    p <- ggplot(aes(x = data[, varX], y = data[, varY], col = group), data = data) + geom_point(size = 1.2, alpha = 0.4)
    if(coefData[varX, "sign"] == "Yes"){
    p <- p + geom_abline(aes(intercept = coefData["Intercept", "mean"], slope = coefData[varX, "mean"]), 
                         lty = 1, lwd = 1.2, data = coefData, color = "black")
    p <- p + xlab(paste(varX)) + ylab(paste(varY))
    }
  } else {
    # Group models
    if(!is.na(str_extract(mdl.name, "G"))){
      # Coefficients
      varXGym <- paste0(varX, ":groupGymnosperms")
      coefData_gr <- mcmcRegDataFun(mod, vars = c(varX, "groupGymnosperms"))
      gimInt <- coefData_gr["Intercept", "mean"] + coefData_gr["groupGymnosperms", "mean"]
      data <- completeFun(pointsData, c(varY, varX))
      # Angiosperms
      p <- ggplot(aes(x = data[, varX], y = data[, varY], col = group), data = data) + geom_point(size = 1.2, alpha = 0.4)
      if(coefData_gr[varX, "sign"] == "Yes"){
      p <- p + geom_abline(aes(intercept = coefData_gr["Intercept", "mean"], slope = coefData_gr[varX, "mean"]), 
                           lty = 1, lwd = 1.2, data = coefData_gr, color = "#e64b35ff")
      }
      p <- p + xlab(paste(varX)) + ylab(paste(varY))
      # Gymnosperm intercept
      if(coefData_gr["groupGymnosperms", "sign"] == "Yes"){
        p <- p + geom_abline(aes(intercept = gimInt,
                                 slope = coefData_gr[varX, "mean"]), data = coefData_gr,
                             color="#4dbbd5ff", size = 1.2)
      } # End only gymnosperm intercept
    }
      # Gymnosperms slope
      if(!is.na(str_extract(mdl.name, "GI"))){
        p <- ggplot(aes(x = data[, varX], y = data[, varY], col = group), data = data) + geom_point(size = 1.2, alpha = 0.4)
        if(coefData_gr[varX, "sign"] == "Yes"){
          p <- p + geom_abline(aes(intercept = coefData_gr["Intercept", "mean"], slope = coefData_gr[varX, "mean"]), 
                               lty = 1, lwd = 1.2, data = coefData_gr, color = "#e64b35ff")
        }
        p <- p + xlab(paste(varX)) + ylab(paste(varY))
        # Coefficients
        varXGym <- paste0(varX, ":groupGymnosperms")
        coefData_gr <- mcmcRegDataFun(mod, vars = c(varX, "groupGymnosperms", varXGym))
        gimInt <- coefData_gr["Intercept", "mean"] + coefData_gr["groupGymnosperms", "mean"]
        gymSlope <- coefData_gr[varX, "mean"] + coefData_gr[varXGym, "mean"]
        gymSlopeLHDP <- coefData_gr[varX, "low95CI"] + coefData_gr[varXGym, "low95CI"]
        gymSlopeHHDP <- coefData_gr[varX, "high95CI"] + coefData_gr[varXGym, "high95CI"]
        gymSlopeSign <- ifelse(gymSlopeLHDP < 0 && gymSlopeHHDP > 0, "No", "Yes" )
        data <- completeFun(pointsData, c(varY, varX))
        if(gymSlopeSign == "Yes"){
          p <- p + geom_abline(aes(intercept = gimInt, 
                                   slope = gymSlope), col = "#4dbbd5ff", lty = 1, lwd = 1, data = coefData_gr)
        }
    }
  } # End model with groups
  p <- p + scale_color_npg() + theme_minimal() + 
    theme(legend.position = "none") + theme( axis.line = element_line(colour = "darkgrey", 
                                                                      size = 1, linetype = "solid"))
  p <- p + theme(axis.title.x =element_blank(), axis.title.y= element_blank()) + xlim(xmin, xmax)
  if(isTRUE(ylim)){ 
    p <- p + ylim(ymin, ymax)
  }
  # Non phylogenetic regression
  noPhyData <-  mcmcRegDataFun(np_mdl, vars = c(varX))
  if(noPhyData[varX, "sign"] == "Yes"){
  p <- p + geom_abline(aes(intercept = noPhyData["Intercept", "mean"], slope = noPhyData[varX, "mean"]),
                       lty = "dashed", lwd = 0.8, data = noPhyData, color = "darkgrey")
  }
  plot(p, res = 1000)
}

bestUniMdlsSumFun <- function(type = "Ht_Env_1"){
  rslts <- list()
  phyInd.res <- data.frame()
  R2.res <- data.frame()
  R2 <- read.csv(paste0(output.dir, "R2_uniModels.csv"), header = T)
  phyInd.df <- read.csv(paste0(output.dir, "phylo_signal_uniModels.csv"), header = T)
  for(i in 1:length(bestUniMdls$type)){
    vars <- unlist(strsplit(as.character(bestUniMdls[i, "variables"]), ", "))[c(1:2)]
    vrs <- paste0(vars[1], "_", vars[2])
    mdl.name <- as.character(paste0(as.character(bestUniMdls[i, "type"]), "_", as.character(bestUniMdls[i, "group"]), "_", vrs))
    smry <- uniMod.res[[mdl.name]]$summary
    smry <- smry$solutions
    rslts[[mdl.name]]$summary <- smry
    write.csv(smry, paste0(manuscript.dir, "figures/Ht_Env/bestUniModels/", mdl.name, ".csv"), row.names = T)
    print(paste0("====>",manuscript.dir, "figures/Ht_Env/bestUniModels/" ,mdl.name, ".csv"))
    # PhyloIndexes best models
    phyInd <- unique(phyInd.df[which(phyInd.df$model == mdl.name), ])
    phyInd.res <- rbind(phyInd.res, phyInd)
    # EVoCor best models
    r2 <- R2[which(R2$model == mdl.name), ]
    R2.res <- rbind(R2.res, r2)
  }
  rslts$r2 <- R2.res
  rslts$phyInd.res <- phyInd.res
  write.csv(R2.res, paste0(manuscript.dir, "figures/Ht_Env/bestUniModels/r2.csv"), row.names = F)
  print(paste0("====>",manuscript.dir, "figures/Ht_Env/bestUniModels/r2.csv"))
  write.csv(phyInd.res, paste0(manuscript.dir, "figures/Ht_Env/bestUniModels/phyloIndexes.csv"), row.names = F)
  print(paste0("====>",manuscript.dir, "figures/Ht_Env/bestUniModels/phyloIndexes.csv"))
  return(rslts)
}


## Phylogenetic reconstruction
genusMapFun <- function(dta = ht.gen, vrble, txn = "genus", tree = ht.tree){
  ## Data
  rec.data <- completeFun(dta, c(paste(vrble), "Genus"))
  rec.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rec.data[, txn]])
  rec.data <- rec.data[rec.data$Genus %in% rec.tree$tip.label, ]
  # Reduction parameter
  nreduction <- 1/solve(length(rec.data$Genus), 60)
  ## Tips
  tips <- rec.tree$tip.label
  tips <- tips[seq(1,length(tips), by= nreduction)]
  ## Vector data
  recVec <- rec.data[, vrble]
  names(recVec) <- rec.data[, txn]
  ## Reconstruction by ML
  contTree <- contMap(tree = rec.tree, x = recVec, plot = F, res = 1000, type = "fan")
  plot(contTree)
  ## get x & y limits for replotting
  lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  xlim <- lastPP$x.lim/7
  ylim <- lastPP$y.lim/7
  # Set colours
  contTree <- setMap(contTree, colorRampPalette(brewer.pal(11, "RdBu"))(diff(range(recVec))))
  #Plot legend
  plot(contTree, ftype = "off", type="fan", lwd=1.8, outline = F)
  # Plot selected tips
  plot(contTree, ftype = "off", type="fan", lwd=1.8, outline = F, legend = F, xlim = xlim, ylim = ylim)
  lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  xpos <- lastPP$xx[sapply(tips,function(x,y) which(y==x),
                           y=rec.tree$tip.label)]
  ypos<-lastPP$yy[sapply(tips,function(x,y) which(y==x),
                         y=rec.tree$tip.label)]
  for(i in 1:length(tips)){
    angle<-atan(ypos[i]/xpos[i])*180/pi
    if(xpos[i]<0) text(xpos[i],ypos[i],paste(gsub("_"," ",tips[i])," ",
                                             sep=""),pos=2,font=3,srt=angle,offset=0,cex=0.8)
    else text(xpos[i],ypos[i],paste(" ",gsub("_"," ",tips[i]),sep=""),
              pos=4,font=3,srt=angle,offset=0,cex= 0.8)
  }
}

famMapFun <- function(dta = ht.gen, vrble, txn = "genus", tree = ht.tree){
  ## Data
  rec.data <- completeFun(dta, c(paste(vrble), txn))
  rec.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rec.data[, txn]])
  rec.data <- rec.data[rec.data[, txn] %in% rec.tree$tip.label, ]
  # Reduction parameter
  # nreduction <- 1/solve(length(rec.data$genus), 60)
  ## Tips
  tips <- rec.tree$tip.label
  # tips <- tips[seq(1,length(tips), by= nreduction)]
  ## Vector data
  recVec <- rec.data[, vrble]
  names(recVec) <- rec.data[, txn]
  ## Reconstruction by ML
  contTree <- contMap(tree = rec.tree, x = recVec, plot = T, res = 1000, type = "fan")
  ## get x & y limits for replotting
  # lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  # xlim <- lastPP$x.lim * 1.5
  # ylim <- lastPP$y.lim * 1.5
  # Set colours
  contTree <- setMap(contTree, colorRampPalette(brewer.pal(11, "Spectral"))(diff(range(recVec))))
  #Plot legend
  plot(contTree, ftype = "off", type="circle", lwd=1.8, outline = F)
  plot(contTree, ftype = "off", type="circle", lwd=1.8, outline = F, legend = F)
  # plot families
  famPlot <- data.frame()
  singlGenPlot <- data.frame()
  for(i in 1:length(rec.data[, "family"])){
    fmly <- rec.data[which(rec.data[, "family"] == rec.data[i, "family"]), c("family", txn)]
    fam.nod <- tips[which(tips %in% fmly$Genus)]
    if(isTRUE(length(fam.nod) == 1)){
      tips.labels <- fam.nod
      singl.Nod.df <- data.frame(tips, tips.labels)
      singlGenPlot <- rbind(singlGenPlot, singl.Nod.df)
    } else {
      nodes <- findMRCA(rec.tree, fam.nod)
      labels <- unique(fmly["family"])
      fam.Nod.df <- data.frame(nodes, labels)
      famPlot <- rbind(famPlot, fam.Nod.df)
    }
  }
  # for(i in 1:length(famPlot$nodes)){
  #   arc.cladelabels(text=as.character(famPlot$family[i]), node=famPlot$nodes[i], 
  #                   orientation = "horizontal", mark.node=F, cex = 0.4)
  # }
  plot(contTree, ftype = "off", type="circle", lwd=1.8, outline = F, legend = F)
  # for(i in 1:length(famPlot$nodes)){
  #   arc.cladelabels(text=NULL, node=famPlot$nodes[i], orientation = "horizontal", mark.node=F, cex = 0.4)
  # }
}

famMapFun <- function(dta = ht.gen, vrble, txn = "Genus", tree = ht.tree){
  ## Data
  rec.data <- completeFun(dta, c(paste(vrble), txn))
  rec.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rec.data[, txn]])
  rec.data <- rec.data[rec.data$Genus %in% rec.tree$tip.label, ]
  # # Reduction parameter
  # nreduction <- 1/solve(length(rec.data$Genus), 60)
  ## Tips
  tips <- rec.tree$tip.label
  # tips <- tips[seq(1,length(tips), by= nreduction)]
  ## Vector data
  recVec <- rec.data[, vrble]
  names(recVec) <- rec.data[, txn]
  ## Reconstruction by ML
  contTree <- contMap(tree = rec.tree, x = recVec, plot = F, res = 1000, type = "fan")
  plot(contTree)
  ## get x & y limits for replotting
  # lastPP <- get("last_plot.phylo",envir=.PlotPhyloEnv)
  # xlim <- lastPP$x.lim * 1.5
  # ylim <- lastPP$y.lim * 1.5
  # Set colours
  contTree <- setMap(contTree, colorRampPalette(brewer.pal(11, "Spectral"))(diff(range(recVec))))
  #Plot legend
  plot(contTree, ftype = "off", type="fan", lwd=1.8, outline = F)
  plot(contTree, ftype = "off", type="fan", lwd=1.8, outline = F, legend = F)
  # plot families
  famPlot <- data.frame()
  singlGenPlot <- data.frame()
  for(i in 1:length(rec.data[, "family"])){
    fmly <- rec.data[which(rec.data[, "family"] == rec.data[i, "family"]), c("family", txn)]
    fam.nod <- tips[which(tips %in% fmly$Genus)]
    if(isTRUE(length(fam.nod) == 1)){
      tips.labels <- fam.nod
      singl.Nod.df <- data.frame(tips, tips.labels)
      singlGenPlot <- rbind(singlGenPlot, singl.Nod.df)
    } else {
      nodes <- findMRCA(rec.tree, fam.nod)
      labels <- unique(fmly["family"])
      fam.Nod.df <- data.frame(nodes, labels)
      famPlot <- rbind(famPlot, fam.Nod.df)
    }
  }
  for(i in 1:length(famPlot$nodes)){
    arc.cladelabels(text=as.character(famPlot$family[i]), node=famPlot$nodes[i],
                    orientation = "horizontal", mark.node=F, cex = 0.4)
  }
  plot(contTree, ftype = "off", type="fan", lwd=1.8, outline = F, legend = F)
  for(i in 1:length(famPlot$nodes)){
    arc.cladelabels(text=NULL, node=famPlot$nodes[i], orientation = "horizontal", mark.node=F, cex = 0.4)
  }
}


phyloMapFun <- function(vars, genus = T, family = T){
  for(var in vars){
  if(isTRUE(genus)){
    f.out <- paste0(manuscript.dir, "figures/reconstructions/genus/gen_phyloMap_", var, ".pdf")
    pdf(f.out, width=10, height=10)
    genusMapFun(ht.gen, var, txn = "Genus", ht.tree)
    paste("======> Image saved on ", f.out)
    dev.off()
  }
  if(isTRUE(family)){
    f.out <- paste0(manuscript.dir, "figures/reconstructions/family/fam_phyloMap_", var, ".pdf")
    pdf(f.out, width=10, height=10)
    famMapFun(ht.gen, var, txn = "Genus", ht.tree)
    paste("======> Image saved on ", f.out)
    dev.off()
  }
  }
}

# Non-group model results
envEffFun <- function(){
  NGmod.res <- list()
  phylo <- read.csv(paste0(output.dir, "phylo_signal_multiModels.csv"), header = T)
  corEvo <- read.csv(paste0(output.dir, "EvoCorr_multiModels.csv"), header = T)
  for(i in 1:length(multiMdls$type)){
    mdl <- multiMod.res[[i]]$name
    if(!is.na(str_extract(mdl, "env_NG"))){
      envEff.df <- as.data.frame(multiMod.res[[i]]$summary$solutions)
      NGmod.res[[mdl]] <- envEff.df
      write.csv(envEff.df, paste0(manuscript.dir, "figures/summary/","envEff_",mdl, ".csv"), row.names = T)
      print(paste0("====>", manuscript.dir, "figures/summary/", "envEff_",mdl, ".csv"))
      envEff.df <- as.data.frame(multiMod.res[[i]]$summary$solutions)
      NGmod.res$envEff <- envEff.df
    }
  }
  evoCor.df <- corEvo[which(!is.na(str_extract(corEvo$model, "NG"))), ]
  write.csv(evoCor.df, paste0(manuscript.dir, "figures/summary/NG_models_evolutionary_correlation.csv"), row.names = F)
  print(paste0("====>", manuscript.dir, "figures/summary/NG_models_evolutionary_correlation.csv"))
  NGmod.res$evoCor.df <- evoCor.df
  phylo.df <- phylo[which(!is.na(str_extract(phylo$model, "NG"))), ]
  write.csv(phylo.df, paste0(manuscript.dir, "figures/summary/NG_models_variance_partition.csv"), row.names = F)
  print(paste0("====>", manuscript.dir, "figures/summary/NG_models_evolutionary_correlation.csv"))
  NGmod.res$phylo <- phylo.df
  return(NGmod.res)
}

# Wittaker 
witPlotFun <- function(wit.dat = ht, vars){
  # Intersection with phylogeny
  phylo <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% wit.dat$genus])
  wit.dat <- wit.dat[wit.dat$genus %in% phylo$tip.label, ]
  # Complete cases
  wit.dat <- completeFun(wit.dat, vars)
  
  wit.dat$MAP_cm <- wit.dat$AP/10
  wit.dat$ID <- rep(1:length(wit.dat$Species))
  my_outliers <- get_outliers(tp = wit.dat[, c("MAT", "MAP_cm")])
  wit.dat$status <- ifelse(wit.dat$ID %in% my_outliers$row_idx, "out", "in")
  print(as.character(wit.dat[which(wit.dat$status == "out"), "Species"]))
  name <- paste(vars, collapse ="_")
  
  p <- whittaker_base_plot() 
  p <- p + geom_point(data = wit.dat, aes(x = MAT, y = MAP_cm, colour = status), size   = 2, shape  = 16, 
                      stroke = 1, alpha  = 0.7) + theme_bw() + theme(legend.background = element_blank())
  leg <- get_legend(p)
  p <- p + theme(legend.position = "none", axis.title =  element_blank()) 
  f.out <- paste0("manuscript/figures/Whittaker/whittaker_", name, ".pdf")
  pdf(f.out, width=3, height=3)
  plot(p, res = 1000)
  paste("======> Image saved on ", f.out)
  dev.off()
  # legend
  f.out <- paste0("manuscript/figures/Whittaker/legend.pdf")
  pdf(f.out, width=5, height=5)
  plot(as_ggplot(leg))
  paste("======> Image saved on ", f.out)
  dev.off()
}

#### RECONSTRUCTIONS (Fig. 1, 2, S5, S6) ----------------------------------------------------------------------------------------------------####

# Data
ht.gen <- read.csv("data/HydraTRY_Genus_2020.csv", header = T)
row.names(ht.gen) <- ht.gen$Genus
ht.tree <- read.tree("data/hydratry_genus.tre") # Genus

phyloMapFun(vars = c("HSM"), genus = F, family = T)


#### Ht ~ Env Plots (Fig. 3) ----------------------------------------------------------------------------------------------------####

uniMdls_DIC <- read.csv(paste0(output.dir, "descriptionAndDIC_uniModels.csv"), header = T)
bestUniMdls <- bestFun(mdls = uniMdls_DIC,  mdlsType = "Ht_Env_1")
write.csv(bestUniMdls, paste0(manuscript.dir, "figures/Ht_Env/bestUniModels/bestUniMdls_Ht_env.csv"), row.names = F)
print(paste0("====> ", manuscript.dir, "figures/Ht_Env/bestUniModels/bestUniMdls_Ht_env.csv"))
bestUniMdlsSum <- bestUniMdlsSumFun()
for(i in 1:length(bestUniMdls$type)){
  # P model
  mdl.name <- paste0(bestUniMdls[i, "type"], "_", bestUniMdls[i, "group"], "_", bestUniMdls[i, "variables"])
  mdl.name <- str_replace(mdl.name, ", ", "_")
  mod <- uniMod.res[[mdl.name]]
  # NP model
  np_mdl.name <- paste0("NP_", mdl.name)
  np_mdl <- NP_uniMod.res[[np_mdl.name]]
  # Variables
  vrs <- as.vector(unlist(strsplit(as.character(bestUniMdls[i, "variables"]), ", ")))
  varY <- vrs[1]
  varX <- vrs[2]
  f.out <- paste0(manuscript.dir, "figures/Ht_Env/", mdl.name, ".pdf")
  pdf(f.out, width= 2, height=2)
  if(!is.na(str_extract(varY, "HSM"))){
    mcmcRegPlotFun(pointsData = plmm.data, xmin = -8, xmax = 8, ylim = T, ymin = -10, ymax = 10)
  } else {
    if(!is.na(str_extract(varY, "Kl"))){
      mcmcRegPlotFun(pointsData = plmm.data, xmin = -8, xmax = 8, ylim = T, ymin = -15, ymax = 0)
    } else {
      mcmcRegPlotFun(pointsData = plmm.data, xmin = -8, xmax = 8, ylim = T, ymin = -4, ymax = 4)
    }
  }
  dev.off()
  cat("==>", f.out, "\n")
}



#### MODELS SUMMARY (Fig. 4) --------------------------------------------------------------------------------------------------####

# Genus level phylogeny
var.part <- read.csv(paste0(output.dir, "phylo_signal_uniModels.csv"), header = T)
write.csv(var.part, paste0(manuscript.dir, "figures/summary/variance_partition.csv"), row.names = F)
print(paste0("====> ", manuscript.dir, "figures/summary/variance_partition.csv"))

# Genus level phylogeny significant evolutionary correlations
signEvoCor.df <- read.csv(paste0(output.dir, "SignEvoCorr_multiModels.csv"), header = T)
write.csv(signEvoCor.df, paste0(manuscript.dir, "figures/summary/Significant_EvolutionaryCorrelations.csv"), row.names = F)
print(paste0("====> ", manuscript.dir, "figures/summary/Significant_EvolutionaryCorrelations.csv"))

# Species level phylogeny significant evolutionary correlations
signEvoCor.df <- read.csv(paste0(output.dir, "extension/Spp_EvoCorr_multiModels.csv"), header = T)
write.csv(signEvoCor.df, paste0(manuscript.dir, "figures/summary/Spp_significant_EvolutionaryCorrelations.csv"), row.names = F)
print(paste0("====> ", manuscript.dir, "figures/summary/Spp_significant_EvolutionaryCorrelations.csv"))

### Environmental effects

envEff <- envEffFun()


#### PRINCIPAL COMPONENTS MAPPING (Fig. S2) -------------------------------------------------------------------------------------####

mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
f.out <- paste0(manuscript.dir, "figures/pca/geo_PC1.pdf")
pdf(f.out, width=36/3, height=18/3)
ggplot(ht) + mapWorld + geom_point(aes(x = longitude_mean, y = latitude_mean, col = PC1))
dev.off()
f.out <- paste0(manuscript.dir, "figures/pca/geo_PC1.pdf")
pdf(f.out, width=36/3, height=18/3)
p <- ggplot(ht) + mapWorld + geom_point(aes(x = longitude_mean, y = latitude_mean, col = PC1))
print(p)
dev.off()
print(paste("==>", f.out))
f.out <- paste0(manuscript.dir, "figures/pca/geo_PC2.pdf")
pdf(f.out, width=36/3, height=18/3)
ggplot(ht) + mapWorld + geom_point(aes(x = longitude_mean, y = latitude_mean, col = PC2))
dev.off()
print(paste("==>", f.out))
f.out <- paste0(manuscript.dir, "figures/pca/geo_PC3.pdf")
pdf(f.out, width=36/3, height=18/3)
ggplot(ht) + mapWorld + geom_point(aes(x = longitude_mean, y = latitude_mean, col = PC3))
dev.off()
print(paste("==>", f.out))

#### WHITTAKER PLOT (Fig. S3) ----------------------------------------------------------------------------------------------------####

witPlotFun(wit.dat = ht, vars = c("PC1", "log_Ks", "log_Hv"))
witPlotFun(wit.dat = ht, vars = c("PC1", "log_Ks", "log_negP50"))
witPlotFun(wit.dat = ht, vars = c("PC1", "log_Ks", "log_negMinWP_md"))
witPlotFun(wit.dat = ht, vars = c("PC1", "log_Hv", "log_negP50"))
witPlotFun(wit.dat = ht, vars = c("PC1",  "log_Hv", "log_negMinWP_md"))
witPlotFun(wit.dat = ht, vars = c("PC1", "log_negMinWP_md", "log_negP50"))
witPlotFun(wit.dat = ht, vars = "PC1")


#### PCA BIPLOT (Fig. S4)---------------------------------------------------------------------------------------------------------####

## Biplot

# PC1 -  PC2
f.out <- paste0(manuscript.dir, "figures/pca/Env_pca_Hydraulic_groups.pdf")
pdf(f.out, width=5, height=5)
pcaPlot(pca, axs = c(1,2))
dev.off()
print(paste("==>", f.out))

# PC1- PC3
f.out <- paste0(manuscript.dir, "figures/pca/Env_pca_Hydraulic_groups_pc3.pdf")
pdf(f.out, width=5, height=5)
pcaPlot(pca, axs = c(1,3))
dev.off()
print(paste("==>", f.out))


## Correlations

# Hydraulic traits PC1-2
f.out <- paste0(manuscript.dir, "figures/pca/Hyd_corr_pc1_pc2.pdf")
pdf(f.out, width=5, height=5)
PcaCorPlot(pca, axs = c(1, 2))
dev.off()
print(paste("==>", f.out))

# Hydraulic traits PC1-3
f.out <- paste0(manuscript.dir, "figures/pca/Hyd_corr_pc1_pc3.pdf")
pdf(f.out, width=5, height=5)
PcaCorPlot(pca, axs = c(1, 3))
dev.off()
print(paste("==>", f.out))


#### GYMNOSPERMS ENVIRONEMNT (Fig. S7)

gym.df <- plmm.data[which(plmm.data$group == "Gymnosperms"), ]
gym.df$fmly <- as.character(gym.df$fmly)

f.out <- paste0(manuscript.dir, "figures/Gymnosperms/Ymin_PC1.pdf")
pdf(f.out, width=3, height=3)
p <- ggplot(gym.df) + geom_point(aes(x = PC1, y = log_negMinWP_md, col = fmly)) + theme(legend.position = "none")
p
dev.off()
print(paste("==>", f.out))

f.out <- paste0(manuscript.dir, "figures/Gymnosperms/HSM_PC1.pdf")
pdf(f.out, width=3, height=3)
p <- ggplot(gym.df) + geom_point(aes(x = PC1, y = HSM, col = fmly)) + theme(legend.position = "none")
p
dev.off()
print(paste("==>", f.out))

f.out <- paste0(manuscript.dir, "figures/Gymnosperms/HSM_PC2.pdf")
pdf(f.out, width=3, height=3)
p <- ggplot(gym.df) + geom_point(aes(x = PC2, y = HSM, col = fmly)) + theme(legend.position = "none")
p
dev.off()
print(paste("==>", f.out))

f.out <- paste0(manuscript.dir, "figures/Gymnosperms/P50_PC1.pdf")
pdf(f.out, width=3, height=3)
p <- ggplot(gym.df) + geom_point(aes(x = PC1, y = log_negP50, col = fmly)) + theme(legend.position = "none")
p
dev.off()
print(paste("==>", f.out))

f.out <- paste0(manuscript.dir, "figures/Gymnosperms/legend.pdf")
pdf(f.out, width=3, height=3)
p <- ggplot(gym.df) + geom_point(aes(x = PC1, y = HSM, col = fmly))
p
dev.off()
print(paste("==>", f.out))
