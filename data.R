library(lipdR)
library(tidyverse)
library(purrr)
library(glue)
library(googlesheets4)
library(geoChronR)


source("functions.R")

set.seed(141)
# Load in conversion ------------------------------------------------------
conv <- read_sheet(ss = "1A7fsNb5SeaUNQMsO8KRZYkubLGEYPJVTJAPaZBOPW1A",sheet = "paleo")
convChron <- read_sheet(ss = "1A7fsNb5SeaUNQMsO8KRZYkubLGEYPJVTJAPaZBOPW1A",sheet = "chron")


# Load in data -------------------------------------------------------


#get reconstructions from all regions
recon_files <- list.files("data/paleo/",pattern = "_reconstruction",full.names = T)
allReg <- map_dfr(recon_files, read_delim)

#manipulate into sites
allSites <- unique(allReg$`ID (Site)`)

#get chronData

s1 <- read_delim("data/chron/Table-S1_chronological_control_points_metadata.csv")
s1$Author_Comment <- str_remove_all(string = s1$Author_Comment,pattern = '[^A-Za-z0-9. ()<>-]')
s1$Material_Dated <- str_remove_all(string = s1$Material_Dated,pattern = '[^A-Za-z0-9. ()<>-]')


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
#sum(!allSites %in% all_cmt_Sites)

#load in other files.
s2 <- read_delim("data/chron/Table-S2_prior_information_of_dates.csv")
s3 <- read_delim("data/chron/Table-S3_bacon_parameter_settings.csv")
#s4 <- read_delim("data/chron/Table-S4_original_chronology_metadata_by_pollen_records.csv")
#s5 <- read_delim("data/chron/Table-S5_legacyage_1_0_chronology.csv")



D <- map2(siteMeta$Site_ID,siteMeta$Dataset_ID,leg2lipd)
good <- which(map_lgl(D,is.list))
Dg <- D[good]
tst <- extractTs(Dg)
geoChronR::plotSummaryTs(tst)


#duplicated?
dn <- map_chr(Dg,"dataSetName")
dup <- duplicated(dn)
dupNames <- unique(dn[dup])

allDup <- which(dn %in% dupNames)

nD <- list()
for(ii in seq_along(dupNames)){
  dni <- which(dn == dupNames[ii])

  nD[[dupNames[ii]]] <- Dg[[dni[1]]]
  for(i in 2:length(dni)){
    nD[[dupNames[ii]]]$paleoData[[i]] <-  Dg[[dni[i]]]$paleoData[[1]]
  }

}


nDts <- extractTs(nD)
dupTs <- extractTs(Dg[allDup])

if(length(nDts) != length(dupTs)){
  stop("these lengths should match")
}

Dgnd <- Dg[-allDup]
an <- map_chr(Dgnd,"dataSetName")
if(length(an) != length(unique(an))){
  stop("should all be unique")
}

names(Dgnd) <- an

Df <- append(Dgnd,nD)


#final check

check <- map_lgl(Df,lipdR:::validLipd)

if(sum(!check) > 0){
  stop("errors")
}

writeLipd(Df,path = "~/Dropbox/lipdverse/LegacyClimate/")





#
# DUP <- Dg[ii]
# ts <- extractTs(DUP) %>% ts2tibble()
# tstemp <- filter(ts,paleoData_primaryTimeseries == TRUE) %>% filter(paleoData_variableName == "temperature")
#
# plotTimeseriesStack(tstemp,time.var = "age")
#
# a <- 1
# plot(tstemp$age[[a]],tstemp$paleoData_values[[a]])
# lines(tstemp$age[[a+2]],tstemp$paleoData_values[[a+2]])

# all_cmt_Sites <- unique(s1$Site_ID)



#oneSite <- filter(allReg,`ID (Site)` ==  allSites[156])
