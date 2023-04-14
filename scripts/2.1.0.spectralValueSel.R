# machine <- "drone"
the_year <- "2019"
data_day <- c("0802", "0810", "0817", "0824", "0831", "0903")

options(stringsAsFactors = FALSE)
source("scripts/1.functionCode.R")
library(readr)
library(stringr)
library(doParallel)

# make the folder
valueFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.0.spectralValueSel"
if (!file.exists(valueFolder)) {
  dir.create(valueFolder)
}

# set the condition
condition <- c("W1", "W2", "W3", "W4")

# choose the VIs name
VIsName <- c("GRVI", "NDVI", "NDRE", "NDI", "RTVI")

cl <- makeCluster(15)
registerDoParallel(cl)
allResult <- foreach(dayInd = 1:length(data_day), .packages = c("doParallel", "stringr")) %dopar% {
  # dayInd <- 1
  the_day <- data_day[dayInd]
  
  spectralEachDayFolder <- paste0("2019_M100_Xacti4eyeCamera_Images/raw/2019", 
                                  the_day)
  refBoardList <- list.files(spectralEachDayFolder, pattern = "ref", full.names = T)
  
  msEachCondition <- foreach(conditionEach = condition) %do% {
    # conditionEach <- condition[2]
    refBoardFile <- refBoardList[grep(conditionEach, refBoardList)]
    refBoardValue <- read.csv(refBoardFile, header = T, row.names = 1)
    colnames(refBoardValue) <- c("475550850_475", "475550850_550", "475550850_850", 
                                 "550660850_550", "550660850_660", "550660850_850", 
                                 "725")
    spectralCsvList0 <- list.files(spectralEachDayFolder, 
                                   full.names = T, pattern = conditionEach)
    spectralCsvList0 <- spectralCsvList0[-grep("ref", spectralCsvList0)]
    # get VIs value
    allMsValue <- lapply(spectralCsvList0, function(eachPlot) {
      # eachPlot <- spectralCsvList0[[1]]
      eachPlotValue <- GetMsData(eachPlot, refBoardValue, VIsName)
      if (ncol(eachPlotValue) < 10) { 
        eachMsData <- matrix(NA, nrow = 4, ncol = nrow(eachPlotValue))
        colnames(eachMsData) <- rownames(eachPlotValue)
      } else {
        eachMsData <- apply(eachPlotValue, 1, function(eachBandValue) {
          # eachBandValue <- eachPlotValue[1, ]
          eachBandMedian <- median(eachBandValue)
          
          maxValue <- max(eachBandValue)
          minValue <- min(eachBandValue)
          interval <- (maxValue - minValue) / round(sqrt(length(eachBandValue)))
          histEach <- hist(eachBandValue, 
                           breaks = seq(minValue, maxValue, interval))
          eachBandMode <- histEach$breaks[which.max(histEach$counts)]
          eachBandQ25 <- quantile(eachBandValue, probs = 0.25)
          eachBandQ75 <- quantile(eachBandValue, probs = 0.75)
          
          return(c(eachBandMedian, eachBandMode, eachBandQ25, eachBandQ75))
        })
      }
      eachPlotMedian <- eachMsData[1, ]
      eachPlotMode <- eachMsData[2, ]
      eachPlotQ25 <- eachMsData[3, ]
      eachPlotQ75 <- eachMsData[4, ]
      
      return(list(eachPlotMedian = eachPlotMedian, eachPlotMode = eachPlotMode, 
                  eachPlotQ25 = eachPlotQ25, eachPlotQ75 = eachPlotQ75))
    })
    
    eachPlotMedianList <- lapply(allMsValue, function(eachMsValue) {
      # eachMsValue <- allMsValue[[1]]
      eachMsValue$eachPlotMedian
    })
    eachPlotModeList <- lapply(allMsValue, function(eachMsValue) {
      # eachMsValue <- allMsValue[[1]]
      eachMsValue$eachPlotMode
    })
    eachPlotQ25List <- lapply(allMsValue, function(eachMsValue) {
      # eachMsValue <- allMsValue[[1]]
      eachMsValue$eachPlotQ25
    })
    eachPlotQ75List <- lapply(allMsValue, function(eachMsValue) {
      # eachMsValue <- allMsValue[[1]]
      eachMsValue$eachPlotQ75
    })
    
    msMedian <- do.call(what = rbind, args = eachPlotMedianList)
    msMode <- do.call(what = rbind, args = eachPlotModeList)
    msQ25 <- do.call(what = rbind, args = eachPlotQ25List)
    msQ75 <- do.call(what = rbind, args = eachPlotQ75List)
    
    plotNumeber <- str_sub(spectralCsvList0, start = 57, 62)
    rownames(msMedian) <- rownames(msMode) <- plotNumeber
    rownames(msQ25) <- rownames(msQ75) <- plotNumeber
    
    return(list(msMedian = msMedian, msMode = msMode, 
                msQ25 = msQ25, msQ75 = msQ75))
  }
  
  medianList <- lapply(msEachCondition, function(msEachConditionList) {
    msEachConditionList$msMedian
  })
  modeList <- lapply(msEachCondition, function(msEachConditionList) {
    msEachConditionList$msMode
  })
  Q25List <- lapply(msEachCondition, function(msEachConditionList) {
    msEachConditionList$msQ25
  })
  Q75List <- lapply(msEachCondition, function(msEachConditionList) {
    msEachConditionList$msQ75
  })
  
  msMedianEachDay <- do.call(what = rbind, args = medianList)
  msModeEachDay <- do.call(what = rbind, args = modeList)
  msQ25EachDay <- do.call(what = rbind, args = Q25List)
  msQ75EachDay <- do.call(what = rbind, args = Q75List)
  
  write.csv(msMedianEachDay, paste0(valueFolder, "/", the_day, "_median.csv"))
  write.csv(msModeEachDay, paste0(valueFolder, "/", the_day, "_mode.csv"))
  write.csv(msQ25EachDay, paste0(valueFolder, "/", the_day, "_Q25.csv"))
  write.csv(msQ75EachDay, paste0(valueFolder, "/", the_day, "_Q75.csv"))
}
stopCluster(cl)


