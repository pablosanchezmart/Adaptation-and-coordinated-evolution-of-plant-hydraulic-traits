################################# TRAITS EVOLUTION RECONSTRUCTION ###################################################################################################


# Pablo Sanchez

### PACKAGES ####
# remove.packages("phytools")
# install.packages("png")

library(phytools)
library(dplyr)
library(pdf)
# library(geiger)
# library(dplyr)
# library(ggtree)

# Represent different clades (families)
#http://blog.phytools.org/2017/03/clade-labels-on-circular-fan-tree.html

# By descendants tips
#http://blog.phytools.org/2017/05/arccladelabels-using-names-of.html

### FUNCTIONS ###

#Function to get complete cases:return subset without NA
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

phyloMapFun <- function(dta, vrble, txn = "Genus", pltTxn, tree){
   ## Data
  rec.data <- completeFun(dta, c(paste(vrble), "Genus"))
  rec.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rec.data[, txn]])
  rec.data <- rec.data[rec.tree$tip.label, ]
  length(rec.data$Genus)
  ## Tips
  tips <- rec.tree$tip.label
  ## Vector data
  recVec <- rec.data[, vrble]
  names(recVec) <- rec.data[, txn]
  ## Reconstruction by ML
  contTree <- contMap(tree = rec.tree, x = recVec, plot = F, res = 1000)
  plot(contTree, ftype = "off", leg.txt = vrble, type="fan", lwd=1.8, outline = F)
  # errorBars <- errorbar.contMap(contTree)
  ## Family plot data (family nodes)
  famPlot <- data.frame()
  singlGenPlot <- data.frame()
  for(i in 1:length(rec.data[, pltTxn])){
    fmly <- rec.data[which(rec.data[, pltTxn] == rec.data[i, pltTxn]), c(pltTxn, "Genus")]
    fam.nod <- tips[which(tips %in% fmly$Genus)]
    if(isTRUE(length(fam.nod) == 1)){
      tips.labels <- fam.nod
      singl.Nod.df <- data.frame(tips, tips.labels)
      singlGenPlot <- rbind(singlGenPlot, singl.Nod.df)
    } else {
      nodes <- findMRCA(rec.tree, fam.nod)
      labels <- unique(fmly[pltTxn])
      fam.Nod.df <- data.frame(nodes, labels)
      famPlot <- rbind(famPlot, fam.Nod.df)
    }
  }
  ## Results
  results <- list(rec.tree, famPlot, contTree, tips)
  return(results)
}


genusMapFun <- function(dta, vrble, txn = "Genus", tree){
  ## Data
  rec.data <- completeFun(dta, c(paste(vrble), "Genus"))
  rec.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rec.data[, txn]])
  rec.data <- rec.data[rec.tree$tip.label, ]
  # Reduction parameter
  nreduction <- 1/solve(length(rec.data$Genus), 60)
  ## Tips
  tips <- rec.tree$tip.label
  ## Vector data
  recVec <- rec.data[, vrble]
  names(recVec) <- rec.data[, txn]
  ## Reconstruction by ML
  contTree <- contMap(tree = rec.tree, x = recVec, plot = F, res = 1000)
  plot(contTree, ftype = "off", type="fan", lwd=1.8, outline = F)
  ## Genus level plot
  tips <- tips[seq(1,length(tips),by= nreduction)]
  ## get x & y limits for replotting
  plot(contTree, type="fan", lwd=1.8, outline = F)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  xlim<-lastPP$x.lim
  ylim<-lastPP$y.lim
  plot(contTree, xlim=xlim, ylim=ylim, ftype = "off", type="fan",lwd=1.8, outline = F)
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


fancyFun <- function(dta, vrble, txn = "Genus", tree){
  ## Data
  rec.data <- completeFun(dta, c(paste(vrble), "Genus"))
  rec.tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% rec.data[, txn]])
  rec.data <- rec.data[rec.tree$tip.label, ]
  ## Vector data
  recVec <- rec.data[, vrble]
  names(recVec) <- rec.data[, txn]
  ## Reconstruction by ML
  fantree <- fancyTree(tree = rec.tree, type="phenogram95", x = recVec)
  return(fantree)
}

# ## Add wBars
# 
# plotTree.wBars(obj$tree,exp(svl),method="plotSimmap",
#                tip.labels=TRUE,fsize=0.7,colors=obj$cols,type="fan",scale=0.002)
# add.color.bar(1.0,obj$cols,title="trait value",lims=obj$lims,prompt=FALSE,
#               x=0.9*par()$usr[1],y=0.9*par()$usr[3])


### DATA ####

setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")

# RData

load(file = "RData/contMap.RData")

ht <- read.csv("Ht_Pablo/HydraTRY_Genus.csv", header = T)
row.names(ht) <- ht$Genus
ht.tree <- read.tree("Phylogeny_data/Ht_tree/hydratry_genus.tre") # Genus

setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis")


### Ks plot ####

## Order plot

ksMap.p <- phyloMapFun(dta = ht, vrble = "log_Ks", txn = "Genus", pltTxn = "order", tree = ht.tree)

for(i in 1:length(ksMap.p[[2]]$nodes)){ 
  arc.cladelabels(text=as.character(ksMap.p[[2]]$order[i]), node=ksMap.p[[2]]$nodes[i], 
                  orientation = "horizontal", mark.node=F, cex = 0.4)
}

## Ks genus

f.out <- paste0("results/log_Ks_phyloMap.pdf")
pdf(f.out, width=10, height=10)
genusMapFun(ht, "log_Ks", txn = "Genus", ht.tree)
paste("======> Image saved on ", f.out)
dev.off()

## Error bars
# errorbar.contMap(ksMap.p)

## Fancytree
# KsFancy.p <- fancyFun(dta = ht, vrble = "log_Ks", txn = "Genus", tree = ht.tree)
# KsFancy.p


### P50 ####
f.out <- paste0("results/log_negP50_phyloMap.pdf")
pdf(f.out, width=10, height=10)
genusMapFun(ht, "log_negP50", txn = "Genus", ht.tree, nreduction = 7)
paste("======> Image saved on ", f.out)
dev.off()

## ERROR BARS
# errorbar.contMap(P50Map.p)

# # KS - P50 PLOT
# rec.data <- completeFun(ht, c("log_Ks", "log_negP50"))
# rec.tree <- drop.tip(ht.tree, ht.tree$tip.label[!ht.tree$tip.label %in% rec.data[, "Genus"]])
# rec.data <- rec.data[rec.tree$tip.label, ]
# data <- as.matrix(rec.data)
# layout(matrix(1:3, 1, 3), widths = c(0.44, 0.12, 0.44))
# par(cex = 1)
# plot(ksMap.p, lwd = 1, outline = F)
# ylim<-c(1-0.12*(length(rec.tree$tip.label)-1),length(rec.tree$tip.label))
# plot.new(); plot.window(xlim=c(-0.1,0.1),ylim=ylim)
# text(rep(0,length(rec.tree$tip.label)), 1:length(rec.tree$tip.label),rec.tree$tip.label)
# plot(P50Map.p, direction="leftwards", sig=1, legend=1)


### log_Hv ####

f.out <- paste0("results/log_Hv_phyloMap.pdf")
pdf(f.out, width=10, height=10)
genusMapFun(ht, "log_Hv", txn = "Genus", ht.tree, nreduction = 9.5)
paste("======> Image saved on ", f.out)
dev.off()

## ERROR BARS
# errorbar.contMap(HvMap.p)

### minWP ####

f.out <- paste0("results/log_negMinWP_md_phyloMap.pdf")
pdf(f.out, width=10, height=10)
genusMapFun(ht, "log_negMinWP_md", txn = "Genus", ht.tree, nreduction = 3.5)
paste("======> Image saved on ", f.out)
dev.off()

## ERROR BARS
# errorbar.contMap(minWPMap.p)

#### SAVE IMAGE ####
save.image(file = "Data/RData/contMap.RData")

