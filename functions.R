
#function to create one file:
leg2lipd <- function(siteId,datasetId){


  # Get Site Metadata -------------------------------------------------------
  thisCmt <- filter(s1,`Site_ID` == siteId & `Dataset_ID` == datasetId)

  cmtMeta <- select(thisCmt,Event:First_Publication) %>%
    distinct()

  if(nrow(cmtMeta) == 1){
    siteMeta <- cmtMeta
  }else if(nrow(cmtMeta) > 1){
    #look at the paleodata to see if they have more lat longs present.
    thisPmt <- filter(allReg,`ID (Site)` == siteId & `ID (Dataset)` == datasetId) %>%
      select(Event:Continent) %>%
      distinct()

    #see if any matches in paleo

    siteMeta <- filter(cmtMeta,near(`Longitude (DD)`,thisPmt$Longitude,tol = 1) & near(`Latitude (DD)`,thisPmt$Latitude,tol = 1) )
    if(nrow(siteMeta) == 0){
      print(glue("no paleo matches for site: {siteId}, dataset: {datasetId}. Skipping."))
      return(NA)
    }else if(nrow(siteMeta) > 1){
      print(glue("multiple paleo matches for site: {siteId}, dataset: {datasetId}. Skipping."))
    }else{
      print(glue("Found 1 paleo match for site: {siteId}, dataset: {datasetId}. Yay?"))

    }
  }else{
    return(NA)
  }


  # Create root -------------------------------------------------------------

  L <- list()
  L$lipdVersion <- 1.3
  L$archiveType <- siteMeta$Archive_Type
  firstAuthor <- str_split(siteMeta$First_Publication,pattern = ",")[[1]][1] %>%
    str_remove_all("[^A-Za-z]")
  if(is.na(firstAuthor)){
    firstAuthor <- "author"
  }
  pubYear <-  str_extract(siteMeta$First_Publication,pattern = "([0-9][0-9][0-9][0-9])")
  if(is.na(pubYear)){
    pubYear <- "1111"
  }

  sitename <- siteMeta$Site_Name %>% str_remove_all(" ") %>%
    str_remove_all("[^A-Za-z0-9]")

  L$dataSetName <- paste(sitename,firstAuthor,pubYear,sep = ".")

  #Be aware there will be duplicate dataset names
  L$datasetId <- paste0("LC",digest::digest(siteMeta))

  L$createdBy <- "github.com/nickmckay/LegacyClimate2LiPD"
  L$originalDataSource <- siteMeta$Data_Source[1]
  L$LegacyClimateSiteId <- siteMeta$Site_ID[1]
  L$LegacyClimateDatasetId <- siteMeta$Dataset_ID[1]

  L <- lipdverseR::initializeChangelog(L)

  # Geo metadata ------------------------------------------------------------
  geo <- list()
  geo$latitude <- siteMeta$`Latitude (DD)`
  geo$longitude <- siteMeta$`Longitude (DD)`
  geo$description <- siteMeta$Site_Description
  geo$siteName <- siteMeta$Site_Name
  geo$continent <- siteMeta$Continent

  L$geo <- geo

  # Pub metadata ------------------------------------------------------------

  pub <- vector(mode = "list",length = 1)
  pub[[1]]$citation <- siteMeta$First_Publication

  L$pub <- pub

  # Paleodata ---------------------------------------------------------------

  #check paleo metadata
  paleoMeta <- filter(allReg,`ID (Site)` == siteId & `ID (Dataset)` == datasetId) %>%
    select(Event:Continent) %>%
    distinct()

  if(nrow(paleoMeta) == 0){
    print("no matches")
    return(NA)
  }

  if(nrow(paleoMeta) > 1){
    stop("this is bad")
  }

  thisPmt <- filter(allReg,`ID (Site)` == siteId & `ID (Dataset)` == datasetId) %>%
    select(contains("["))

  an <- names(thisPmt)

  pmt <-vector(mode = "list",length = 1)

  for(coll in 1:length(an)){
    #get conversion row
    cr <- which(conv$ColName == an[coll])
    if(length(cr) != 1){
      stop("bad conv match")
    }
    thisCol <- list()
    thisCol$number <- coll
    thisCol$TSid <- lipdR::createTSid("lc")
    thisCol$values <- as.matrix(thisPmt[,coll])
    thisCol$variableName <- conv$variableName[cr]
    thisCol$units <- conv$units[cr]

    if(startsWith(thisCol$variableName,"age")){
      thisCol$values <- as.matrix(thisPmt[,coll])*1000
      thisCol$units <- "yr BP"
    }

    thisCol$primaryTimeseries <- conv$primaryTimeseries[cr]
    thisCol$primaryAgeColumn <- conv$primaryAgeColumn[cr]
    thisCol$summaryStatistic <- conv$`summary statistic`[cr]
    if(!is.na(conv$calibration_method[cr])){
      thisCol$calibration <- list()
      thisCol$calibration$method <- conv$calibration_method[cr]
      thisCol$calibration$target <- conv$calibration_target[cr]
      thisCol$interpretation <- vector(mode = "list",length = 1)
      thisCol$interpretation[[1]]$scope <- "climate"
      thisCol$interpretation[[1]]$variable <- conv$interpVariable[cr]
      thisCol$interpretation[[1]]$variableDetail <- conv$interpVariableDetail[cr]
      thisCol$interpretation[[1]]$seasonality <- conv$seasonality[cr]
    }

    cname <- conv$variableName[cr]
    ni <- 0
    while(any(names(pmt[[1]]) == cname)){
      ni <- ni+1
      cname <- paste0(conv$variableName[cr],"-",ni)
    }
    pmt[[1]][[cname]] <- thisCol
  }

  L$paleoData <- vector(mode = "list",length = 1)
  L$paleoData[[1]]$measurementTable <- pmt

  # chronData ---------------------------------------------------------------


  #check chron metadata
  chronMeta <- thisCmt %>%
    select(Event:Continent) %>%
    distinct()

  if(nrow(chronMeta) != 1){
    thisCmt <- filter(thisCmt,near(`Longitude (DD)`,geo$longitude,tol = 1) & near(`Latitude (DD)`,geo$latitude,tol = 1))

    chronMeta <- thisCmt %>%
      select(Event:Continent) %>%
      distinct()
    if(nrow(chronMeta) != 1){
      print("this is bad for chron")
      return(NA)
    }
  }

  #prep data frame
  elev <- unique(thisCmt$`Elevation (m a.s.l.)`)
  if(length(elev) == 1){
    if(!is.na(elev)){
      L$geo$elevation <- elev
    }
  }

  cd <- thisCmt %>%
    select(!Event:Continent) %>%
    mutate(age14C = `Age_Uncalibrated (kyr BP)`*1000) %>%
    mutate(age = `Age_Calibrated (kyr BP)`*1000) %>%
    mutate(ageUncertainty = `Calibrated dating_Error (kyr)`*1000)

  cd$age14CUncertainty <- 1000*(cd$`Dating Error_Older (kyr)` + cd$`Dating Error_Younger (kyr)`)/2
  cd$depth_top <- cd$`Depth (cm)` - cd$`Thickness (cm)`/2
  cd$depth_bottom <- cd$`Depth (cm)` + cd$`Thickness (cm)`/2

  toRem <- c("Age_Uncalibrated (kyr BP)", "Dating Error_Older (kyr)","Dating Error_Younger (kyr)","Age_Calibrated (kyr BP)","Calibrated dating_Error (kyr)")
  cd <- cd %>%
    select(-c(1:4)) %>%
    select(!all_of(toRem))


  cdn <- names(cd)

  cmt <-vector(mode = "list",length = 1)

  for(i in 1:length(cdn)){
    #get conversion row
    cr <- which(convChron$ColName == cdn[i])
    if(length(cr) != 1){
      stop("bad conv match")
    }
    thisCol <- list()
    thisCol$number <- i
    thisCol$TSid <- lipdR::createTSid("lc-ch")
    thisCol$values <- as.matrix(cd[,i])
    thisCol$variableName <- convChron$variableName[cr]
    thisCol$units <- convChron$units[cr]
    if(!is.na(convChron$`summary statistic`[cr])){
      thisCol$summaryStatistic <- convChron$`summary statistic`[cr]
    }


    cmt[[1]][[thisCol$variableName]] <- thisCol
  }

  L$chronData <- vector(mode = "list",length = 1)
  L$chronData[[1]]$measurementTable <- cmt


  return(L)
}
