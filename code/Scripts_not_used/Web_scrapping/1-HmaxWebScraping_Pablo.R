####################### Routines to automate web scraping from various sites ###
####################### M Mencuccini August 2018 ###

## libraries (some may not be needed) ####
#install.packages("phyndr", dependencies = T)
library(grDevices)
library(Matrix)
library(lattice)
library(graphics)
library(aptg) # aptg: Automatic Phylogenetic Tree Generator to test if spp is in angiosperm tree; takes too long!!!
library(Taxonstand) # step 1: clean up nomenclature
library(taxonlookup) # step 2: look it up  <-- problems with the 3.5.1 R vesion
library(phyndr) # step 3: swap species in tree for others with traits <-- problems with the 3.5.1 R vesion
library(plyr)
library(phytools)
library(ape)
library(caper)
library(dplyr)
library(stringr)
library(data.table)
library(rlang)
library(xml2)
library(rvest)
library(httr)
library(selectr)

setwd("C:/Users/pablo/Desktop/Master/TFM/Analisis/Data/Ht_Pablo")
ht1 <- read.csv("HydraTRY_Pablo_original.csv", header = T)

## Hmax webscrapping
# Function for numbers extraction
numextract <- function(string){ 
  str_extract_all(string, "\\d+\\.*\\d*")                    # I've change the pattern because we'll don't have negative nubmers. Original pattern: "\\-*\\d+\\.*\\d*"
} 

# Function for Hmax extraction in the loop
Hmax.extract <- function(x, p){
  if(any(grepl(p, x))){
    str <- str_extract_all(x, p)                                  # Sentence containing Hmax
    Hmax <- numextract(str)                           # Hmax extraction using numextract (previously defined)
    return(Hmax)
  } else {                                                        # No Hmax record for the specie
    (NA)
  } # end if
} # end function


## Example
text <- c(" sdf (15-60 cm  in height) with a spreading habit.")
numextract(text)
pattern <- "-\\d+\\.*\\d*"
Hmax.extract(text, pattern)


### SCRAPPING OF MISSING Hmax ###

# SCRAPING EFLORA of CHINA for Hmax data ####
ht <- read.csv("HydraTRY_Pablo_HmaxCompl.csv", header = T)

missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

#Checking new species
# ht[ht$Species == "Oryza sativa", "Hmax"]

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  session <- html_session("http://www.efloras.org")
  form <- html_form(session)[[1]]
  form <- set_values(form, name_str=NameStr)                # fills the form on main page
  new_url <- submit_form(session,form)                      # gets session details with address of new page
  
  new_session <- html_session(new_url$url)                  # starts new session using spp. name
  text <- new_session %>% html_nodes("#TableMain") %>% html_text
  
  if(is.na(str_extract(text, "No taxa found"))){
    n <- n + 1
    present.sp[n] <- paste(NameStr)
    cat(paste(NameStr), "\n")
  }
}

length(present.sp) # 161 spp

## Patterns (ordered from more to less reliable)

# Metres
p1 <- "to \\d+\\.*\\d*\\ m tall"                        # Check 28 new
p2 <- "\\-*\\d+\\.*\\d*\\ m tall"                       # Check 30 new
p3 <- "\\-*\\(*\\d+\\.*\\d*\\)*\\ m "                   # Check, 6 new
p4 <- "\\(\\-\\d+\\)\\ m tall"                          # Check, 0 new
p5 <- "\\(\\-\\d+\\)\\ m"                               # Check, 3 new
p6 <- "to \\d+\\.*\\d*\\ m"                             # Check, 1 new (carefull with this one, not only height)
p7 <- "\\[\\-*\\d+\\.*\\d*\\]\\ m tall"                 
p8 <- "\\-\\d+\\.*\\d*\\ m"                             # Check, 2 new

#dm
p1dm <- "to \\d+\\.*\\d*\\ dm tall"                       # Check, 0 new
p2dm <- "\\-*\\d+\\.*\\d*\\ dm tall"                      # Check, 0 new
p3dm <- "\\(\\-\\d+\\)\\ dm tall" 
p4dm <- " \\d+\\.*\\d*\\ dm." 

#cm
p1cm <- "to \\d+\\.*\\d*\\ cm tall"                       # Check, 0 new
p2cm <- "\\-*\\d+\\.*\\d*\\ cm tall"                      # Check, 0 new
p3cm <- "\\-*\\d+\\.*\\d*\\ cm"                           #Check, 0 new
p4cm <- "\\-*\\d+\\ cm"                                   # Check, 2 new

pat.vec <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9)
pat.vec_dm <- c(p1dm, p2dm, p3dm, p4dm)
pat.vec_cm <- c(p1cm, p2cm, p3cm, p4cm)

## Extraction bucle
Hmax.values <- vector()
several.values <- character()
n <- 0

#metres
for(p in pat.vec) {                                        # Bucle for scrapping all the patterns
cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    session <- html_session("http://www.efloras.org")
    form <- html_form(session)[[1]]
    form <- set_values(form, name_str=NameStr)                # fills the form on main page
    new_url <- submit_form(session,form)                      # gets session details with address of new page
    
    new_session <- html_session(new_url$url)                  # starts new session using spp. name
    
    # saves body of message and extracts record
    text <- new_session %>% html_nodes("body") %>% html_text
    
    # processes string
    numbers = numextract(text)[[1]]
    sp_numb = numbers[nchar(numbers)>3][1]                    # eliminate some shortish vectors and looks at the first one
    
    #correct final url
    # new_url = paste("http://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=", sp_numb, sep = "")
    new_url = paste("http://www.efloras.org/florataxon.aspx?flora_id=1&taxon_id=", sp_numb, sep = "")   # Sometimes it changes to this
    
    # new session at correct location in Flora of China  
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(paste(NameStr), "\n")
        n + 1
        several.values[n] <- paste(NameStr) 
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])
      } else {
        if(is.na(Hmax.values[i])){
          Hmax.values[i] <- NA 
        }
      } # End of if extraction
  } # End of scrapping bucle
cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# decimeters
for(p in pat.vec_dm) {                                        # Bucle for scrapping all the patterns 
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    session <- html_session("http://www.efloras.org")
    form <- html_form(session)[[1]]
    form <- set_values(form, name_str=NameStr)                # fills the form on main page
    new_url <- submit_form(session,form)                      # gets session details with address of new page
    
    new_session <- html_session(new_url$url)                  # starts new session using spp. name
    
    # saves body of message and extracts record
    text <- new_session %>% html_nodes("body") %>% html_text
    
    # processes string
    numbers = numextract(text)[[1]]
    sp_numb = numbers[nchar(numbers)>3][1]                    # eliminate some shortish vectors and looks at the first one
    
    #correct final url
    # new_url = paste("http://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=", sp_numb, sep = "")
    new_url = paste("http://www.efloras.org/florataxon.aspx?flora_id=1&taxon_id=", sp_numb, sep = "")   # Sometimes it changes to this
    
    # new session at correct location in Flora of China  
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(paste(NameStr), "\n")
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])/10
    } else {
      if(is.na(Hmax.values[i])){
        Hmax.values[i] <- NA 
      }  
    } # End of if extraction
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# centimeters
for(p in pat.vec_cm) {                                        # Bucle for scrapping all the patterns 
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    session <- html_session("http://www.efloras.org")
    form <- html_form(session)[[1]]
    form <- set_values(form, name_str=NameStr)                # fills the form on main page
    new_url <- submit_form(session,form)                      # gets session details with address of new page
    
    new_session <- html_session(new_url$url)                  # starts new session using spp. name
    
    # saves body of message and extracts record
    text <- new_session %>% html_nodes("body") %>% html_text
    
    # processes string
    numbers = numextract(text)[[1]]
    sp_numb = numbers[nchar(numbers)>3][1]                    # eliminate some shortish vectors and looks at the first one
    
    #correct final url
    # new_url = paste("http://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id=", sp_numb, sep = "")
    new_url = paste("http://www.efloras.org/florataxon.aspx?flora_id=1&taxon_id=", sp_numb, sep = "")   # Sometimes it changes to this
    
    # new session at correct location in Flora of China  
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(paste(NameStr), "\n")
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])/100
    } else {
      if(is.na(Hmax.values[i])){
        Hmax.values[i] <- NA 
      } 
    } # End of if extraction
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# Bind species and values
length(Hmax.values[which(!is.na(Hmax.values))])                 # to check how many were recovered
missing.spp <- as.data.frame(missing)
colnames(missing.spp) <- "missing.spp" 
scraped_Hmax <- as.data.frame(cbind(missing.spp, Hmax.values))
scraped_Hmax <- scraped_Hmax[which(!is.na(scraped_Hmax$Hmax.values)), ]
scraped_Hmax

## Add values to HydraTRY dataset

new.Hmax <- data.frame(1)
#new.Hmax <- read.csv("NewHmaxValues.csv", header = T)
n <- 0  # Counter, we'll use it later

for(i in 1:length(scraped_Hmax$missing.spp)){                                                                                    # Bucle for all the missing spp
  sp <- scraped_Hmax$missing.spp[i]                                                                                              # Select a specie of the missing spp set
  if(is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){    # HydraTRY: NA, new.Hmax: no NA (Hmax value)
  ht[which(ht$Species == sp), ]$Hmax <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values                        # Assign new Hmax value to HtdraTRY
  cat("Hmax for", scraped_Hmax$missing.spp[i], "scrapped.\n")                                                                      # Print that we have a new record
  } else {                                                                                                                       # If HydraTRY: value and/or new.Hmax: NA
    if(!is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){ # If we have a new value for a specie that already had Hmax in HydraTRY
      if(ht[which(ht$Species == sp), ]$Hmax != scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values) {              # If Hmax values are different
      cat("New measure for", scraped_Hmax$missing.spp[i], ".\nHmax Ht =", ht[which(ht$Species == sp), ]$Hmax, ".\n Scrapped =",  # Print the new measure
         scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values, ".\n")                               
      n <- n + 1                                                                                                                 # Counter for the record of the new measures
      new.Hmax[n, "Species"] <- sp
      new.Hmax[n, "Hmax"] <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values
    } # end new measure
    } else {                                                                                                                     # If both values are the same, next.
    next()                                                      
    }
  } # End repeated measure
} # End of bucle

new.Hmax

missing.sppn <- with(ht, Species[which(is.na(Hmax))])
length(missing.sppn)

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)


# SCRAPING PLANTS FOR A FUTURE ####
ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "+")
  #correct final url
  new_url = paste("https://pfaf.org/user/plant.aspx?LatinName=", NameStr, sep = "")
  # new session 
  new_session <- html_session(new_url)                          # starts new session using spp. name
  
  text <- new_session %>% html_nodes("#ctl00_ContentPlaceHolder1_lbldisplatinname") %>% html_text
  NameStr <- as.character(missing[i])
  if(!is.na(str_extract(text, NameStr))){
    n <- n + 1
    present.sp[n] <- NameStr
    cat(NameStr, "\n")
  }
}

length(present.sp)

## Patterns (ordered from more to less reliable)

# meters
p1 <- "to \\d+\\.*\\d*\\ m tall"         # Check, 3 new
p2 <- "\\-*\\d+\\.*\\d*\\ m tall"        # Ceck, 0 new
p3 <- "\\-*\\(*\\d+\\.*\\d*\\)*\\ m tall"# Ceck, 3 new
p4 <- "growing to \\d+\\.*\\d*\\ m "     # Ceck, 24 new
p5 <- "up to \\d+\\.*\\d*\\ m"             # Ceck, 0 new
p6 <- "up to \\d+\\.*\\d*\\ m "           # Ceck, 0 new
p7 <- "about \\d+\\.*\\d*\\ m "           # Ceck, 0 new

# dm
p1dm <- "to \\d+\\.*\\d*\\ dm tall"        # Ceck, 0 new
p2dm <- "\\-*\\d+\\.*\\d*\\ dm tall"       # Ceck, 0 new

# cm
p1cm <- "to \\d+\\.*\\d*\\ cm tall"        # Ceck, 0 new
p2cm <- "\\-*\\d+\\.*\\d*\\ cm tall"       # Ceck, 0 new

pat.vec <- c(p1, p2, p3, p4, p5, p5, p7, p8)
pat.vec_dm <- c(p1dm, p2dm)
pat.vec_cm <- c(p1cm, p2cm)
Hmax.values <- vector()

## Extraction bucle

# meters
for(p in pat.vec) {                                        # Bucle for scrapping all the patterns 
  cat("scrapping pattern",p)
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "+")
    #correct final url
    new_url = paste("https://pfaf.org/user/plant.aspx?LatinName=", NameStr, sep = "")
    # new session at correct location
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])
    } else {
      if(is.na(Hmax.values[i])){
        Hmax.values[i] <- NA 
      }
    } # End of if extraction
  } # End of scrapping bucle
cat("\nScrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# decimeters
for(p in pat.vec_dm) {                                        # Bucle for scrapping all the patterns 
  cat("scrapping pattern",p)
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "+")
    #correct final url
    new_url = paste("https://pfaf.org/user/plant.aspx?LatinName=", NameStr, sep = "")
    # new session at correct location
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])/10
    } else {
      if(is.na(Hmax.values[i])){
        Hmax.values[i] <- NA 
      }
    } # End of if extraction
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# centimeters
for(p in pat.vec_cm) {                                        # Bucle for scrapping all the patterns 
  cat("scrapping pattern",p)
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "+")
    #correct final url
    new_url = paste("https://pfaf.org/user/plant.aspx?LatinName=", NameStr, sep = "")
    # new session at correct location
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])/100
    } else {
      if(is.na(Hmax.values[i])){
        Hmax.values[i] <- NA 
      }
    } # End of if extraction
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# Bind species and values
length(Hmax.values[which(!is.na(Hmax.values))])                 # to check how many were recovered
missing.spp <- as.data.frame(missing)
colnames(missing.spp) <- "missing.spp" 
scraped_Hmax <- as.data.frame(cbind(missing.spp, Hmax.values))
scraped_Hmax <- scraped_Hmax[which(!is.na(scraped_Hmax$Hmax.values)), ]

## Add values to HydraTRY dataset

new.Hmax <- data.frame(1)
#new.Hmax <- read.csv("NewHmaxValues.csv", header = T)
n <- 0  # Counter, we'll use it later

for(i in 1:length(scraped_Hmax$missing.spp)){                                                                                    # Bucle for all the missing spp
  sp <- scraped_Hmax$missing.spp[i]                                                                                              # Select a specie of the missing spp set
  if(is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){    # HydraTRY: NA, new.Hmax: no NA (Hmax value)
    ht[which(ht$Species == sp), ]$Hmax <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values                        # Assign new Hmax value to HtdraTRY
    cat("Hmax for", scraped_Hmax$missing.spp[i], "scrapped.")                                                                      # Print that we have a new record
  } else {                                                                                                                       # If HydraTRY: value and/or new.Hmax: NA
    if(!is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){ # If we have a new value for a specie that already had Hmax in HydraTRY
      if(ht[which(ht$Species == sp), ]$Hmax != scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values) {              # If Hmax values are different
        cat("New measure for", scraped_Hmax$missing.spp[i], ".\nHmax Ht =", ht[which(ht$Species == sp), ]$Hmax, ".\n Scrapped =",  # Print the new measure
            scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values, ".\n")                               
        n <- n + 1                                                                                                                 # Counter for the record of the new measures
        new.Hmax[n, "Species"] <- sp
        new.Hmax[n, "Hmax"] <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values
      } # end new measure
    } else {                                                                                                                     # If both values are the same, next.
      next()                                                      
    }
  } # End repeated measure
} # End of bucle

new.Hmax

missing.sppn <- with(ht, Species[which(is.na(Hmax))])
length(missing.sppn)

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)


# SCRAPING USEFUL TROPICAL PLANTS ####
ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "+")
  #correct final url
  new_url = paste("http://tropical.theferns.info/viewtropical.php?id=", NameStr, sep = "")
  # new session 
  new_session <- html_session(new_url)                          # starts new session using spp. name
  
  text <- new_session %>% html_nodes(".PageBox:nth-child(1)") %>% html_text
  if(is.na(str_extract(text, "no record"))){
    n <- n + 1
    present.sp[n] <- NameStr
    cat(NameStr, "\n")
  }
}

length(present.sp)

## Patterns (ordered from more to less reliable)

#metres
p1 <- "grow up to \\-*\\d+\\.*\\d*\\ metres tall"     # Check, 25 new
p2 <- "specimens up to \\d+\\.*\\d*\\ metres"         # Check, 2 new
p3 <- "specimens to \\d+\\.*\\d*\\ metres"            # Check, 3 new
p4 <- "exceptionally to \\d+\\.*\\d*\\ metres"        # Check, 2 new 
p5 <- "grow to \\-*\\d+\\.*\\d*\\ metres tall"        # Check,0 new
p6 <- "occasionally to \\d+\\.*\\d*\\ metres"         # Check, 2 new
p7 <- "up to \\d+\\.*\\d*\\ metres tall"              # Check, 5 sp new
p8 <- "\\-*\\ \\d+\\.*\\d*\\ metres tall"             # Check, 51 sp new
p9 <- "\\d+\\.*\\d*\\ metres tall"                    # Check 0 new
p10 <- "\\-*\\ \\d+\\.*\\d*\\ metres"                 # Carefull, just to Duguetia furfuracea (the others are long, not heigh)
p11 <- "reach a height of up to \\d+\\.*\\d*\\ metres"# Check, 1 sp new

#dm

#cm
p1cm <- "growing up to \\d+cm\\ tall"                   # Check, few new
p2cm <- "growing about \\d+cm\\ tall"                  # Check, few new

## Extraction bucle

pat.vec <- c(p1, p2, p3, p4, p5, p5, p7, p8, p9, p10, p11)
pat.vec_cm <- c(p1cm, p2cm)
Hmax.values <- vector()

# metres
for(p in pat.vec) {                                        # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "+")
    #correct final url
    new_url = paste("http://tropical.theferns.info/viewtropical.php?id=", NameStr, sep = "")
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes(".PageBox:nth-child(1)") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
        }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])
    } else {
      Hmax.values[i] <- NA
    } # End of if extraction
  } # End of scrapping bucle
cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# dm

# cm
for(p in pat.vec_cm) {                                        # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "+")
    #correct final url
    new_url = paste("http://tropical.theferns.info/viewtropical.php?id=", NameStr, sep = "")
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes(".PageBox:nth-child(1)") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])/100
    } else {
        if(is.na(Hmax.values[i]))
      Hmax.values[i] <- NA
    } # End of if extraction
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle


# Bind species and values
length(Hmax.values[which(!is.na(Hmax.values))])                 # to check how many were recovered
missing.spp <- as.data.frame(missing)
colnames(missing.spp) <- "missing.spp" 
scraped_Hmax <- as.data.frame(cbind(missing.spp, Hmax.values))
scraped_Hmax <- scraped_Hmax[which(!is.na(scraped_Hmax$Hmax.values)), ]

## Add values to HydraTRY dataset

new.Hmax <- data.frame(1)
#new.Hmax <- read.csv("NewHmaxValues.csv", header = T)
n <- 0  # Counter, we'll use it later

for(i in 1:length(scraped_Hmax$missing.spp)){                                                                                    # Bucle for all the missing spp
  sp <- scraped_Hmax$missing.spp[i]                                                                                              # Select a specie of the missing spp set
  if(is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){    # HydraTRY: NA, new.Hmax: no NA (Hmax value)
    ht[which(ht$Species == sp), ]$Hmax <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values                        # Assign new Hmax value to HtdraTRY
    cat("Hmax for", scraped_Hmax$missing.spp[i], "scrapped.")                                                                      # Print that we have a new record
  } else {                                                                                                                       # If HydraTRY: value and/or new.Hmax: NA
    if(!is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){ # If we have a new value for a specie that already had Hmax in HydraTRY
      if(ht[which(ht$Species == sp), ]$Hmax != scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values) {              # If Hmax values are different
        cat("New measure for", scraped_Hmax$missing.spp[i], ".\nHmax Ht =", ht[which(ht$Species == sp), ]$Hmax, ".\n Scrapped =",  # Print the new measure
            scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values, ".\n")                               
        n <- n + 1                                                                                                                 # Counter for the record of the new measures
        new.Hmax[n, "Species"] <- sp
        new.Hmax[n, "Hmax"] <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values
      } # end new measure
    } else {                                                                                                                     # If both values are the same, next.
      next()                                                      
    }
  } # End repeated measure
} # End of bucle


new.Hmax

missing.sppn <- with(ht, Species[which(is.na(Hmax))])
length(missing.sppn)

# Mannual scrapping

ht[ht$Species == "Pyrrosia adnascens", ]$Hmax <- 0.12

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)


# SCRAPING PLANT NET ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "~")
  #correct final url
  new_url = paste("http://plantnet.rbgsyd.nsw.gov.au/cgi-bin/NSWfl.pl?page=nswfl&lvl=sp&name=", NameStr, sep = "")
  # new session 
  new_session <- html_session(new_url)                          # starts new session using spp. name
  
  text <- new_session %>% html_nodes(".hide+ table , .normal > tbody td , h3") %>% html_text
  if(is.na(str_extract(text, "No entry"))){
    n <- n + 1
    present.sp[n] <- NameStr
    cat(NameStr, "\n")
  }
}

length(present.sp)

## Patterns (ordered from more to less reliable)

#metres
p1 <- "\\-*\\d+\\.*\\d*\\ m high"
p2 <-  "to \\d+\\.*\\d*\\ m high"

## Extraction bucle

pat.vec <- c(p1, p2)
Hmax.values <- vector()

# metres
for(p in pat.vec) {                                        # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "~")
    #correct final url
    new_url = paste("http://plantnet.rbgsyd.nsw.gov.au/cgi-bin/NSWfl.pl?page=nswfl&lvl=sp&name=", NameStr, sep = "")
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes("body") %>% html_text
    
    if(!is.na(Hmax.extract(text, p))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])
    } else {
      Hmax.values[i] <- NA
    } # End of if extraction
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

# Bind species and values
length(Hmax.values[which(!is.na(Hmax.values))])                 # to check how many were recovered
missing.spp <- as.data.frame(missing)
colnames(missing.spp) <- "missing.spp" 
scraped_Hmax <- as.data.frame(cbind(missing.spp, Hmax.values))
scraped_Hmax <- scraped_Hmax[which(!is.na(scraped_Hmax$Hmax.values)), ]

## Add values to HydraTRY dataset

new.Hmax <- data.frame(1)
#new.Hmax <- read.csv("NewHmaxValues.csv", header = T)
n <- 0  # Counter, we'll use it later

for(i in 1:length(scraped_Hmax$missing.spp)){                                                                                    # Bucle for all the missing spp
  sp <- scraped_Hmax$missing.spp[i]                                                                                              # Select a specie of the missing spp set
  if(is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){    # HydraTRY: NA, new.Hmax: no NA (Hmax value)
    ht[which(ht$Species == sp), ]$Hmax <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values                        # Assign new Hmax value to HtdraTRY
    cat("Hmax for", scraped_Hmax$missing.spp[i], "scrapped.")                                                                      # Print that we have a new record
  } else {                                                                                                                       # If HydraTRY: value and/or new.Hmax: NA
    if(!is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){ # If we have a new value for a specie that already had Hmax in HydraTRY
      if(ht[which(ht$Species == sp), ]$Hmax != scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values) {              # If Hmax values are different
        cat("New measure for", scraped_Hmax$missing.spp[i], ".\nHmax Ht =", ht[which(ht$Species == sp), ]$Hmax, ".\n Scrapped =",  # Print the new measure
            scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values, ".\n")                               
        n <- n + 1                                                                                                                 # Counter for the record of the new measures
        new.Hmax[n, "Species"] <- sp
        new.Hmax[n, "Hmax"] <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values
      } # end new measure
    } else {                                                                                                                     # If both values are the same, next.
      next()                                                      
    }
  } # End repeated measure
} # End of bucle

new.Hmax

missing.sppn <- with(ht, Species[which(is.na(Hmax))])
length(missing.sppn)

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)


# SCRAPING WIKIPEDIA ENGLISH ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "_")
  #correct final url
  new_url = paste("https://en.wikipedia.org/wiki/", NameStr, sep = "")
  
  if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
    next()
  } else { # hence wiki page exists
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    text <- new_session %>% html_nodes("body") %>% html_text
    if(is.na(str_extract(text, "does not have an article"))){
      n <- n + 1
      present.sp[n] <- str_replace_all(NameStr, "_", " ")
      cat(NameStr, "\n")
    }
  }
} # End of presence searching

length(present.sp)

## Patterns (ordered from more to less reliable)

#metres
p1 <-  "to \\d+\\.*\\d*\\ m"
p2 <- "-\\d+\\.*\\d*\\ m tall"
p3 <- "about \\d+\\.*\\d*\\ metres"
p4 <- "\\d+\\.*\\d*\\ m"

# cm
p1cm <-  "-\\d+\\.*\\d*"

## Extraction bucle

pat.vec <- c(p1)
pat.vec_cm <- c(p1cm)
Hmax.values <- vector()

# metres
for(p in pat.vec) {                                                # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "_")
    #correct final url
    new_url = paste("https://en.wikipedia.org/wiki/", NameStr, sep = "")
    if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
      Hmax.values[i] <- NA
      next()
    } else { # hence wiki page exists
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    
    text <- new_session %>% html_nodes(".biota+ p") %>% html_text
    
    if(!is.na(Hmax.extract(text, p4))) {                          # this loop makes no scraping; just uses info from prior session
      Hmax.char <- Hmax.extract(text, p4)
      if(str_length(Hmax.char) > 2){ 
        cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
      }
      Hmax.values[i] <- as.numeric(Hmax.char[[1]])
    } else {
      Hmax.values[i] <- NA
  } # End of if extraction
 }
} # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

#cm
for(p in pat.vec_cm) {                                                # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "_")
    #correct final url
    new_url = paste("https://en.wikipedia.org/wiki/", NameStr, sep = "")
    if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
      Hmax.values[i] <- NA
      next()
    } else { # hence wiki page exists
      # new session 
      new_session <- html_session(new_url)                          # starts new session using spp. name
      
      text <- new_session %>% html_nodes(".biota+ p") %>% html_text
      
      if(!is.na(Hmax.extract(text, p1cm))) {                          # this loop makes no scraping; just uses info from prior session
        Hmax.char <- Hmax.extract(text, p1cm)
        if(str_length(Hmax.char) > 2){ 
          cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
        }
        Hmax.values[i] <- as.numeric(Hmax.char[[1]])/100
      } else {
        Hmax.values[i] <- NA
      } # End of if extraction
    }
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle


# Bind species and values
length(Hmax.values[which(!is.na(Hmax.values))])                 # to check how many were recovered
missing.spp <- as.data.frame(missing)
colnames(missing.spp) <- "missing.spp" 
scraped_Hmax <- as.data.frame(cbind(missing.spp, Hmax.values))
scraped_Hmax <- scraped_Hmax[which(!is.na(scraped_Hmax$Hmax.values)), ]


## Add values to HydraTRY dataset

new.Hmax <- data.frame(1)
n <- 0  # Counter, we'll use it later

for(i in 1:length(scraped_Hmax$missing.spp)){                                                                                    # Bucle for all the missing spp
  sp <- scraped_Hmax$missing.spp[i]                                                                                              # Select a specie of the missing spp set
  if(is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){    # HydraTRY: NA, new.Hmax: no NA (Hmax value)
    ht[which(ht$Species == sp), ]$Hmax <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values                        # Assign new Hmax value to HtdraTRY
    cat("Hmax for", scraped_Hmax$missing.spp[i], "scrapped.")                                                                      # Print that we have a new record
  } else {                                                                                                                       # If HydraTRY: value and/or new.Hmax: NA
    if(!is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){ # If we have a new value for a specie that already had Hmax in HydraTRY
      if(ht[which(ht$Species == sp), ]$Hmax != scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values) {              # If Hmax values are different
        cat("New measure for", scraped_Hmax$missing.spp[i], ".\nHmax Ht =", ht[which(ht$Species == sp), ]$Hmax, ".\n Scrapped =",  # Print the new measure
            scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values, ".\n")                               
        n <- n + 1                                                                                                                 # Counter for the record of the new measures
        new.Hmax[n, "Species"] <- sp
        new.Hmax[n, "Hmax"] <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values
      } # end new measure
    } else {                                                                                                                     # If both values are the same, next.
      next()                                                      
    }
  } # End repeated measure
} # End of bucle



new.Hmax

missing.sppn <- with(ht, Species[which(is.na(Hmax))])
length(missing.sppn)

# Manual assignation
# ht[ht$Species == "Alnus alnobetula", ]$Hmax <- 12
# ht[ht$Species == "Bauhinia galpinii", ]$Hmax <- 6
# ht[ht$Species == "Eriogonum cinereum", ]$Hmax <- 1.22
# ht[ht$Species == "Garrya ovata", ]$Hmax <- 1.8
# ht[ht$Species == "Hyptis emoryi", ]$Hmax <- 2.44
# ht[ht$Species == "Lomatia polymorpha", ]$Hmax <- 4
# ht[ht$Species == "Lupinus arizonicus", ]$Hmax <- 0.50
# ht[ht$Species == "Malacothamnus fasciculatus", ]$Hmax <- 5
# ht[ht$Species == "Melanobatus leucodermis", ]$Hmax <- 3
# ht[ht$Species == "Nestegis sandwicensis", ]$Hmax <- 20
# ht[ht$Species == "Orites diversifolia", ]$Hmax <- 3
# ht[ht$Species == "Quercus oleoides", ]$Hmax <- 15
# ht[ht$Species == "Salvia leucophylla", ]$Hmax <- 1.5
# ht[ht$Species == "Sambucus cerulea", ]$Hmax <- 9

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)





# SCRAPING WIKIPEDIA SPANISH ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "_")
  #correct final url
  new_url = paste("https://es.wikipedia.org/wiki/", NameStr, sep = "")
  
  if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
    next()
  } else { # hence wiki page exists
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    text <- new_session %>% html_nodes("body") %>% html_text
    if(is.na(str_extract(text, "no tiene"))){
      n <- n + 1
      present.sp[n] <- str_replace_all(NameStr, "_", " ")
      cat(NameStr, "\n")
    }
  }
} # End of presence searching

length(present.sp)

## Patterns (ordered from more to less reliable)

#metres
p1 <-  "to \\d+\\.*\\d*\\ m"
p2 <- "-\\d+\\.*\\d*\\ m tall"
p3 <- "about \\d+\\.*\\d*\\ metres"
p4 <- "\\d+\\.*\\d*\\ m"

# cm
p1cm <-  "-\\d+\\.*\\d*"

## Extraction bucle

pat.vec <- c(p1)
pat.vec_cm <- c(p1cm)
Hmax.values <- vector()

# metres
for(p in pat.vec) {                                                # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "_")
    #correct final url
    new_url = paste("https://en.wikipedia.org/wiki/", NameStr, sep = "")
    if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
      Hmax.values[i] <- NA
      next()
    } else { # hence wiki page exists
      # new session 
      new_session <- html_session(new_url)                          # starts new session using spp. name
      
      text <- new_session %>% html_nodes(".biota+ p") %>% html_text
      
      if(!is.na(Hmax.extract(text, p4))) {                          # this loop makes no scraping; just uses info from prior session
        Hmax.char <- Hmax.extract(text, p4)
        if(str_length(Hmax.char) > 2){ 
          cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
        }
        Hmax.values[i] <- as.numeric(Hmax.char[[1]])
      } else {
        Hmax.values[i] <- NA
      } # End of if extraction
    }
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle

#cm
for(p in pat.vec_cm) {                                                # Bucle for scrapping all the patterns
  cat("scrapping pattern", p, "\n")
  for (i in 1:length(missing)) {
    uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
    NameStr <- missing[i]
    NameStr <- str_replace_all(NameStr, " ", "_")
    #correct final url
    new_url = paste("https://en.wikipedia.org/wiki/", NameStr, sep = "")
    if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
      Hmax.values[i] <- NA
      next()
    } else { # hence wiki page exists
      # new session 
      new_session <- html_session(new_url)                          # starts new session using spp. name
      
      text <- new_session %>% html_nodes(".biota+ p") %>% html_text
      
      if(!is.na(Hmax.extract(text, p1cm))) {                          # this loop makes no scraping; just uses info from prior session
        Hmax.char <- Hmax.extract(text, p1cm)
        if(str_length(Hmax.char) > 2){ 
          cat(NameStr, "\n")                                         # In case that there are more than one occurrences (check individual cases)
        }
        Hmax.values[i] <- as.numeric(Hmax.char[[1]])/100
      } else {
        Hmax.values[i] <- NA
      } # End of if extraction
    }
  } # End of scrapping bucle
  cat("Scrapped values: ", length(Hmax.values[which(!is.na(Hmax.values))]), ".\n")
} # End of pattern bucle


# Bind species and values
length(Hmax.values[which(!is.na(Hmax.values))])                 # to check how many were recovered
missing.spp <- as.data.frame(missing)
colnames(missing.spp) <- "missing.spp" 
scraped_Hmax <- as.data.frame(cbind(missing.spp, Hmax.values))
scraped_Hmax <- scraped_Hmax[which(!is.na(scraped_Hmax$Hmax.values)), ]


## Add values to HydraTRY dataset

new.Hmax <- data.frame(1)
n <- 0  # Counter, we'll use it later

for(i in 1:length(scraped_Hmax$missing.spp)){                                                                                    # Bucle for all the missing spp
  sp <- scraped_Hmax$missing.spp[i]                                                                                              # Select a specie of the missing spp set
  if(is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){    # HydraTRY: NA, new.Hmax: no NA (Hmax value)
    ht[which(ht$Species == sp), ]$Hmax <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values                        # Assign new Hmax value to HtdraTRY
    cat("Hmax for", scraped_Hmax$missing.spp[i], "scrapped.")                                                                      # Print that we have a new record
  } else {                                                                                                                       # If HydraTRY: value and/or new.Hmax: NA
    if(!is.na(ht[which(ht$Species == sp), ]$Hmax) && !is.na(scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values)){ # If we have a new value for a specie that already had Hmax in HydraTRY
      if(ht[which(ht$Species == sp), ]$Hmax != scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values) {              # If Hmax values are different
        cat("New measure for", scraped_Hmax$missing.spp[i], ".\nHmax Ht =", ht[which(ht$Species == sp), ]$Hmax, ".\n Scrapped =",  # Print the new measure
            scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values, ".\n")                               
        n <- n + 1                                                                                                                 # Counter for the record of the new measures
        new.Hmax[n, "Species"] <- sp
        new.Hmax[n, "Hmax"] <- scraped_Hmax[which(scraped_Hmax$missing.spp == sp), ]$Hmax.values
      } # end new measure
    } else {                                                                                                                     # If both values are the same, next.
      next()                                                      
    }
  } # End repeated measure
} # End of bucle



new.Hmax

missing.sppn <- with(ht, Species[which(is.na(Hmax))])
length(missing.sppn)

# Manual assignation
# ht[ht$Species == "Alberta magna", ]$Hmax <- 13
# ht[ht$Species == "Ipomoea wolcottiana", ]$Hmax <- 12
# ht[ht$Species == "Pereskia guamacho", ]$Hmax <- 8
# ht[ht$Species == "Pereskia marcanoi", ]$Hmax <- 6
# ht[ht$Species == "Pereskia weberiana", ]$Hmax <- 3
# ht[ht$Species == "Retama sphaerocarpa", ]$Hmax <- 3
# ht[ht$Species == "Retama sphaerocarpa", ]$Hmax <- 3

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)





# SCRAPING WIKIPEDIA FRENCH ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "_")
  #correct final url
  new_url = paste("https://fr.wikipedia.org/wiki/", NameStr, sep = "")
  
  if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
    next()
  } else { # hence wiki page exists
    # new session 
      n <- n + 1
      present.sp[n] <- str_replace_all(NameStr, "_", " ")
      cat(NameStr, "\n")
    }
} # End of presence searching

length(present.sp)

# Manual assignation
# ht[ht$Species == "Elymus spicatus", ]$Hmax <- 1
# ht[ht$Species == "Justicia californica", ]$Hmax <- 1.5



write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)

# SCRAPING OAKS OF THE WORLD ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "_")
  NameStr <- tolower(NameStr)
  #correct final url
  new_url = paste("http://oaks.of.the.world.free.fr/", NameStr, ".htm", sep = "")
  # new session 
  if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
    next()
  } else { # hence wiki page exists
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
      n <- n + 1
      present.sp[n] <- NameStr
      cat(NameStr, "\n")
    }
  } # End of presence searching


length(present.sp)  # 2 spp

# ht[ht$Species == "Quercus floribunda", ]$Hmax <- 30
# ht[ht$Species == "Quercus sebifera", ]$Hmax <- 2

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)





# SCRAPING FLORA BASE WESTERN AUSTRALIA ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "+")
  #correct final url
  new_url = paste("https://florabase.dpaw.wa.gov.au/search/quick?q=", NameStr, sep = "")
  # new session 
  new_session <- html_session(new_url)                          # starts new session using spp. name
  
  text <- new_session %>% html_nodes("body") %>% html_text
  if(is.na(str_extract(text, "Nothing was found"))){
    n <- n + 1
    present.sp[n] <- NameStr
    cat(NameStr, "\n")
  }
}

length(present.sp)

## Patterns (ordered from more to less reliable)

# ht[ht$Species == "Christella dentata", ]$Hmax <- 0.8
# ht[ht$Species == "Actinostrobus pyramidalis", ]$Hmax <- 4 

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)

# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)

# SCRAPING AFRICAN PLANTS ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  session <- html_session("https://www.prota4u.org/database/")
  form <- html_form(session)[[1]]
  form <- set_values(form, species = NameStr)                # fills the form on main page
  new_url <- submit_form(session,form)                      # gets session details with address of new page
  
  new_session <- html_session(new_url$url)                  # starts new session using spp. name
  text <- new_session %>% html_nodes("body") %>% html_text
  
  if(is.na(str_extract(text, "No taxa found"))){
    n <- n + 1
    present.sp[n] <- paste(NameStr)
    cat(paste(NameStr), "\n")
  }
}

length(present.sp) 


## Patterns (ordered from more to less reliable)

# ht[ht$Species == "Christella dentata", ]$Hmax <- 0.8

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)

# If there are new Hmax values
# write.csv(new.Hmax, "Web_scrapping/NewHmaxValues.csv",  row.names = FALSE)

# GLOBAL PLANTS (JSTOR) ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0
for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "+")
  #correct final url
  new_url = paste("https://plants.jstor.org/search?filter=name&so=ps_group_by_genus_species+asc&Query=", NameStr, sep = "")
  # new session 
  new_session <- html_session(new_url)                          # starts new session using spp. name
  
  text <- new_session %>% html_nodes("body") %>% html_text
  if(is.na(str_extract(text, "No entry"))){
    n <- n + 1
    present.sp[n] <- NameStr
    cat(NameStr, "\n")
  }
}

length(present.sp)



# Manual assignation
# ht[ht$Species == "Acosmium cardenasii", ]$Hmax <- 4       #https://plants.jstor.org/stable/10.5555/al.ap.specimen.k000500775?searchUri=filter%3Dname%26so%3Dps_group_by_genus_species%2Basc%26Query%3DAcosmium%2Bcardenasii
# ht[ht$Species == "Aporosa globifera", ]$Hmax <- 9.14
write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)

# SCRAPING CASA BIO ####

ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing) 

## Searching missing species 
present.sp <- character()
n <- 0

for (i in 1:length(missing)) {
  uastring <- "Mozilla/5.0 (Windows NT 6.1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2228.0 Safari/537.36"
  NameStr <- missing[i]
  NameStr <- str_replace_all(NameStr, " ", "-")
  NameStr <- tolower(NameStr)
  #correct final url
  new_url = paste("https://casabio.org/taxa/", NameStr, sep = "")
  if(http_error(new_url)) { # avoids error if page is missing (normally Wikipedia has solution)
    next()
  } else { # hence wiki page exists
    # new session 
    new_session <- html_session(new_url)                          # starts new session using spp. name
    n <- n + 1
    present.sp[n] <- NameStr
    cat(NameStr, "\n")
  }
} # End of presence searching

length(present.sp)  # 13 spp

# Manually addition
# ht[ht$Species == "Aspalathus hirta", ]$Hmax <- 2
# ht[ht$Species == "Aspalathus pachyloba", ]$Hmax <- 2
# ht[ht$Species == "Cliffortia ruscifolia", ]$Hmax <- 1.5
# ht[ht$Species == "Erica plukenetii", ]$Hmax <- 1.5
# ht[ht$Species == "Erica plukenetii", ]$Hmax <- 1.5

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)

# SCRAPING OTHER WEBS MANUALLY ####
# write.csv(missing, "last_missing.csv", row.names = F)
ht <- read.csv("HydraTRY_Pablo_logData_HmaxCompl.csv", header = T)
missing <- with(ht, Species[which(is.na(Hmax))])
length(missing)

missing <- read.csv("last_missing.csv", header = T)
missing

# ht[ht$Species == "Actinostemon klotzschii", ]$Hmax <- 10    #http://flora.ipe.org.br/sp/17
# ht[ht$Species == "Adiantum petiolatum", ]$Hmax <- 0.4       #https://biogeodb.stri.si.edu/bioinformatics/croat/specie/Adiantum%20petiolatum,e,n
# ht[ht$Species == "Aglaia glabrata", ]$Hmax <- 20            #http://portal.cybertaxonomy.org/flora-malesiana/cdm_dataportal/taxon/b5b545ad-bf1b-46d6-92ca-a9fc6917ef5f
# ht[ht$Species == "Anacardium nanum", ]$Hmax <- 1.5          #http://sweetgum.nybg.org/science/world-flora/details.php?irn=28158
# ht[ht$Species == "Anacardium nanum", ]$Hmax <- 1.5  
# ht[ht$Species == "Ficus benguetensis", ]$Hmax <- 15           #http://portal.cybertaxonomy.org/flora-malesiana/cdm_dataportal/taxon/3727357d-2c63-4823-afec-2b3836bcbf98
# ht[ht$Species == "Vatica bella", ]$Hmax <- 50                 #http://portal.cybertaxonomy.org/flora-malesiana/cdm_dataportal/taxon/692837dc-11e3-4db5-9376-ca107b8dd93a
# ht[ht$Species == "Salsola iberica", ]$Hmax <- 0.6               #https://www.asturnatura.com/especie/salsola-kali.html

write.csv(ht, "HydraTRY_Pablo_logData_HmaxCompl.csv", row.names = FALSE)
