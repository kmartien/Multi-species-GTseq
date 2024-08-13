library(tidyverse)

samps <- read.csv("data-raw/Ggri_archive_dump_5June2024.csv")
tissue.data <- read.csv("data-raw/Ggri_Tissue_Vial_data.csv")
dna.data <- read.csv("data-raw/Ggri_DNA_Concentration.csv")

names(tissue.data)[1] <- "LABID"
names(dna.data)[1] <- "LABID"

samp.dat <- select(samps, c(LABID, CONTACT_NAME, 
                      Restrict_Description, Date_Collected, Latitude, Longitude, 
                      Location, Island, State, Country, Ocean_Basin))

samp.dat$DQuant <- do.call(rbind, lapply(samp.dat$LABID, function(i){
  d.dat <- filter(dna.data, LABID == i)
  if(length(which(!is.na(d.dat$DQuant))) == 0) return("") else return(paste(d.dat$DQuant, collapse = ","))
}))

samp.dat$concentration <- do.call(rbind, lapply(samp.dat$LABID, function(i){
  d.dat <- filter(dna.data, LABID == i)
  if(length(which(!is.na(d.dat$Concentration_ng_ul))) == 0) return("") else return(paste(d.dat$Concentration_ng_ul, collapse = ", "))
}))

samp.dat$dna.comments <- do.call(rbind, lapply(samp.dat$LABID, function(i){
  d.dat <- filter(dna.data, LABID == i)
  if(length(which(!is.na(d.dat$dbo_vw_DNA_Storage_Comments))) == 0) return("") else return(paste(d.dat$dbo_vw_DNA_Storage_Comments, collapse = ", "))
}))

samp.dat$TQuant <- do.call(rbind, lapply(samp.dat$LABID, function(i){
  d.dat <- filter(tissue.data, LABID == i)
  if(length(which(!is.na(d.dat$TQuant))) == 0) return("") else return(paste(d.dat$TQuant, collapse = ","))
}))

write.csv(samp.dat, file = "qry_results/Ggri.lcWGS.sample.picking.csv")
