############################################### VARIABLE PREPARATION ##############################################

# Pablo Sanchez Martínez
# 18/03/2019


### PACKAGES ####

#install.packages("humidity", dependencies = T)

#vpd 
library(humidity) # for saturation vapor preassure calculation from T
#https://cran.r-project.org/web/packages/humidity/humidity.pdf


setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data")
ht <- read.csv("Ht_Pablo/HydraTRY_Pablo.csv", header = T)

### FUNCTION ####

# To divide temperature variables per 10
tempVarCorrection.fun <- function(df, vars){
  for(i in vars){
    df[[i]] <- df[[i]]/10
  } # end bucle
  return(df)
} # end function

### VARIABLE CHECKING ####
names(ht)

### Change variable names ####
names(ht)[names(ht) == "ALT"] <- "altitude"
names(ht)[names(ht) == "Growth.form"] <- "growth.form"
names(ht)[names(ht) == "Life.form"] <- "life.form"
names(ht)[names(ht) == "Leaf.form"] <- "leaf.form"
names(ht)[names(ht) == "longitude"] <- "longitude_mean"
names(ht)[names(ht) == "latitude"] <- "latitude_mean"
names(ht)[names(ht) == "TMean.year"] <- "MAT"
names(ht)[names(ht) == "PH_HOX_200cm_1km_ll"] <- "ph.200cm"
names(ht)[names(ht) == "ClayContent_percentage_60cm_1km"] <- "clay.60cm_perc"
names(ht)[names(ht) == "SiltContent_percentage_60cm_1km"] <- "silt.60cm_perc"
names(ht)[names(ht) == "SandContent_percentage_60cm_1km"] <- "sand.60cm_perc"
names(ht)[names(ht) == "SaturatedWaterContent_volumetricFraction_200cm_1km"] <- "SWC.200cm"
names(ht)[names(ht) == "AbsoluteDepthBedRock_1km_ll"] <- "ADBRock"
names(ht)[names(ht) == "Prec.year"] <- "MAP"

### Redundant or useless variables ####
ht$MAP_mean <- NULL # MAP
ht$MAT_mean <- NULL # MAP
ht$PredictedClassTAXNWRB_1km_ll <- NULL # since we have the percentage contents.
ht$ai.values <- NULL # AI
ht$PH_HOX_15cm_1kmll <- NULL #ph
ht$long_mean <- NULL # longitude
ht$lat_mean <- NULL # latitude

### VARIABLE CALCULATION/RECALCULATION ####

### vars/10 ####

# For temperature and PH vars we need to divide the values by 10
divVars <- c("MAT", "Diurnal.T.Range", "TMax.warmest", "TMin.Coldest", "TRange",
              "TMean.Wettest", "TMean.Driest", "TMean.Warmest", "TMean.Coldest", "ph.200cm")

ht <- tempVarCorrection.fun(ht, divVars)

## Isothermality/100
ht$Isothermality <- ht$Isothermality/100

## Temperature seasonality/100
ht$Seasonality <- ht$Seasonality/100

## AI*0.0001
ht$AI <- ht$AI * 0.0001

## Some of the variables need to be calculated:

### VPDMax calculation ####

ht$VPDMax <- NA
for(i in 1:length(ht$Species)){
  if(all(is.na(ht[i, 80:91]))){
    ht$VPDMax[i] <- NA
  } else {
  tmax.month <- as.character(which.max(ht[i , 80:91]))                                       # Warmest month
  vapr.month <- paste0("vapr_", tmax.month)                                                  # vapr (kPa) of the warmest month           
  tmax.month <- paste0("tmax_", tmax.month)                                                  # Tmax value warmest month
  actVapr <- ht[i, vapr.month]                                                               # vapr value warmest month

  satVp <- SVP(ht[i, tmax.month], isK = F, formula = c("Clausius-Clapeyron", "Murray"))      # SatVP (hPa) (Tmax warmest month)
  ht$VPDMax[i] <- (satVp/10) - actVapr                                                       # VPDmax = vapr(warmest month) - Tmax (warmest month)
  }
} # end of function

### wind.Max ####
ht$windMax <- NA
for(i in 1:length(ht$Species)){
  ht$windMax[i] <- max(ht[i , 104:115])
}

### srad.min ####
ht$sradMin <- NA
for(i in 1:length(ht$Species)){
  ht$sradMin[i] <- min(ht[i , 92:103])
}

### Delete non-usefull climate data ###
names(ht)
ht <- ht[, -c(55:126)]
names(ht)
summary(ht)

# Save FINAL database
ht <- ht[, ]
"Species", "genus", "family", "order", "group", "growth.form", "life.form", "leaf.form", "biome", "Hmax", "Hact", "Ks", "Kl", "P50", "MinWP_pd", "MinWP_md", "AlAs", "WD", 
"Vdia", "VD", "VLmax", "altitude", "Hv", "HSM", "PItlp", "WUE", "SLA", "Aarea", "Amass", "LL"              
[31] "pLL"              "N"                "Narea"            "Rd"               "longitude_mean"   "latitude_mean"   
[37] "tmin"             "tmax"             "srad"             "wind"             "vapr"             "MAP"             
[43] "MAT"              "Diurnal.T.Range"  "Isothermality"    "Seasonality"      "TMax.warmest"     "TMin.Coldest"    
[49] "TRange"           "TMean.Wettest"    "TMean.Driest"     "TMean.Warmest"    "TMean.Coldest"    "Prec.Wettest"    
[55] "Prec.Driest"      "Prec.Seasonality" "Prec.WetQ"        "Prec.DryQ"        "Prec.WarmQ"       "Prec.ColdQ"      
[61] "ph.200cm"         "clay.60cm_perc"   "silt.60cm_perc"   "sand.60cm_perc"   "SWC.200cm"        "ADBRock"         
[67] "AI"               "VPDMax"           "windMax"          "sradMin"         
> 

write.csv(ht, "HT_Pablo_2019.csv", row.names = F)
