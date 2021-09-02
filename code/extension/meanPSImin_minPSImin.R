#### RELATIONSHIP BETWEEN MEAN AND MIN PSImin ------------------------------------------------------------------------------ ####

library(dplyr)
library(plyr)

# Mean function
mean.fun <- function(x){
  y <- mean(x, na.rm = T)
  if(is.infinite(y)){
    y <- NA
  }
  return(y)
}

# HydraTRY
ht <- read.csv("data/HydraTRY_Pablo_2020.csv", header = T)  # Includes three principal components

# Dissagregated dataset
xf <- readxl::read_excel("data/xylem_functional_traits_database_2014.xls")
xf$Species <- paste0(xf$Genus, " ", xf$Species)
xf$Species <- as.factor(xf$Species)

# Sample size

minWP_N <- xf %>% filter(!is.na(xf$MinWP_md))

length(minWP_N$Species)
length(unique(minWP_N$Species))

tb <- table(minWP_N$Species)

barplot(tb)

length(which(tb == 1))
length(which(tb > 1))
length(which(tb > 5))
length(which(tb > 10))

# Mean PSImin
mean_MinWP_md <- aggregate(xf$MinWP_md, by = list(xf$Species), FUN = mean.fun) 
names(mean_MinWP_md) <- c("Species", "mean_MinWP_md")

# Combine
comb.df <- join(ht[, c("Species", "MinWP_md")], mean_MinWP_md, by = "Species", type = "full")
compl.df <- comb.df[complete.cases(comb.df$MinWP_md, comb.df$mean_MinWP_md), ]
length(compl.df$Species)

# Model
cor <- cor(log(-compl.df$mean_MinWP_md), log(-compl.df$MinWP_md))
cor

mod <- lm(log(-compl.df$mean_MinWP_md) ~ log(-compl.df$MinWP_md))
summary(mod)

plot(x = log(-compl.df$mean_MinWP_md), y = log(-compl.df$MinWP_md))
