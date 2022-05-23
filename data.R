library(lipdR)
library(tidyverse)
library(purrr)
library(glue)



# Load in data -------------------------------------------------------


#get reconstructions from all regions
recon_files <- list.files("data/paleo/",pattern = "_reconstruction",full.names = T)
allReg <- map_dfr(recon_files, read_delim)

#manipulate into sites
allSites <- unique(allReg$`ID (Site)`)

#get chronData

s1 <- read_delim("data/chron/Table-S1_chronological_control_points_metadata.csv")

#apparently you need both site and dataset ids to get unique sites?
siteMeta <- s1 %>%
  group_by(Site_ID,Dataset_ID) %>%
  summarize(nSiteName = length(unique(Site_Name)),.groups = "keep")

#yay!
if(any(sum(siteMeta$nSiteName > 1))){
  stop("not uniquely identifying sites.")
}

#get site/dataset combos


#see if any sites missing chron data
sum(!allSites %in% all_cmt_Sites)

#load in other files.
s2 <- read_delim("data/chron/Table-S2_prior_information_of_dates.csv")
s3 <- read_delim("data/chron/Table-S3_bacon_parameter_settings.csv")
#s4 <- read_delim("data/chron/Table-S4_original_chronology_metadata_by_pollen_records.csv")
#s5 <- read_delim("data/chron/Table-S5_legacyage_1_0_chronology.csv")



#function to create one file:
leg2lipd <- function(siteId,datasetId){
  thisCmt <- filter(s1,`Site_ID` == siteId & `Dataset_ID` == datasetId)

  cmtMeta <- select(thisCmt,Event:First_Publication) %>%
    distinct()



  if(nrow(cmtMeta) > 1){
    #look at the paleodata to see if they have more lat longs present.

    print(glue("I don't think there should be more than 1 row, but there is in site {siteId}"))
  }

  #thisPaleo <- filter(allReg,`ID (Site)` == siteId)



}

walk2(siteMeta$Site_ID,siteMeta$Dataset_ID,leg2lipd)



all_cmt_Sites <- unique(s1$Site_ID)



oneSite <- filter(allReg,`ID (Site)` ==  allSites[156])
