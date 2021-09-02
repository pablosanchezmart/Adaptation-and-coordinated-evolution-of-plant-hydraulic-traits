############################################# TABLES AND FIGURES ########################################################################

#### PACKAGES ####
# install.packages("wesanderson")
library(ggbiplot)
library(wesanderson)
names(wes_palettes)

setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")

#### R DATA ###

## PCA
load("RData/PCA.RData")

## PPCA


## PHYLOGENETIC RECONSTRUCTIONS


## MCMCglmm
load("RData/Phylogenetic_linear_mixed_multi_models.RData")

## PATH ANALYSES

#### FIGURE PCA ####

## Biplot
f.out <- paste0("C:/Users/pablo/Desktop/Master/TFM/Documents/Manuscript/Figures/Env_pca_Hydraulic_groups.png")
png(f.out, width=1610, height=925, units = 'px')
p <- fviz_pca_biplot(pca, axes = c(1,2), label = "none", 
                     geom.ind = "point", 
                     habillage = ht_pca_gr, alpha.ind = 0.5,
                  # addEllipses=TRUE, ellipse.level=0.95, 
                  select.var = list(contrib = 10),
                  col.var = "#3182bd",
                  repel = T,
                  title = "Environemntal PCA - Hydraulic traits",
                  ggtheme = theme(axis.text=element_text(size=10), axis.title=element_text(size=16, face="bold"))
)
p <- p + xlim(-8, 8) + ylim (-8, 8)
p <- p + scale_color_manual(values = c("#99d8c9", "#2ca25f", "#e5f5e0")) +  theme_minimal()
# Add hydraulic traits correlations
p <- fviz_add(p, hyd.coord*10, addlabel = F,
              geom="arrow", repel = T, linetype = "solid", color = "black",
              ggtheme = theme(axis.text=element_text(size=10), axis.title=element_text(size=16, face="bold"))
)
print(p)
dev.off()
print(paste("==>", f.out))
print(p)


