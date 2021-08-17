########## 2.1.0.spectralValueSel ############
GetMsData <- function(plotCsv, refBoardValue, VIsName) {
  plotMat0 <- as.matrix(read.csv(plotCsv, header = T, row.names = 1))
  spectral475_475 <- plotMat0["475550850_475", ] / refBoardValue[, "475550850_475"]
  spectral475_550 <- plotMat0["475550850_550", ] / refBoardValue[, "475550850_550"]
  spectral475_850 <- plotMat0["475550850_850", ] / refBoardValue[, "475550850_850"]
  
  spectral550_550 <- plotMat0["550660850_550", ] / refBoardValue[, "550660850_550"]
  spectral550_660 <- plotMat0["550660850_660", ] / refBoardValue[, "550660850_660"]
  spectral550_850 <- plotMat0["550660850_850", ] / refBoardValue[, "550660850_850"]
  spectral725 <- plotMat0["725", ] / refBoardValue[, "725"]
  
  spectral850 <- (spectral475_850 + spectral550_850) / 2
  # add VIs
  GRVI <- (spectral550_550 - spectral550_660) / (spectral550_550 + spectral550_660)
  NDVI <- (spectral550_850 - spectral550_660) / (spectral550_850 + spectral550_660)
  GNDVI <- (spectral550_850 - spectral550_550) / (spectral550_850 + spectral550_550)
  BNDVI <- (spectral475_850 - spectral475_475) / (spectral475_850 + spectral475_475)
  NDRE <- (spectral850 - spectral725) / (spectral850 + spectral725)
  CIgreen <- (spectral550_850 / spectral550_550) - 1
  CIRE <- (spectral850 / spectral725) - 1
  MSR <- ((spectral550_850 / spectral550_660) - 1) / sqrt((spectral550_850 / spectral550_660) + 1)
  MSRRE <- ((spectral850 / spectral725) - 1) / sqrt((spectral850 / spectral725) + 1)
  NDI <- (spectral725 - spectral550_660) / (spectral725 + spectral550_660)
  BGVI <- (spectral475_550 - spectral475_475) / (spectral475_550 + spectral475_475)
  EVI <- 2.5 * (spectral850 - spectral550_660) / (spectral850 + 6 * spectral550_660 - 7.5 * spectral475_475 + 1)
  VARI <- (spectral550_550 - spectral550_660) / (spectral550_550 + spectral550_660 - spectral475_475)
  ARI <- (1 / spectral550_550) - (1 / spectral725)
  mARI <- ((1 / spectral550_550) - (1 / spectral725)) * spectral550_850
  
  BRVI <- (spectral475_475 - spectral550_660) / (spectral475_475 + spectral550_660)
  BRENDI <- (spectral725 - spectral475_475) / (spectral725 + spectral475_475)
  GRENDI <- (spectral725 - spectral550_550) / (spectral725 + spectral550_550)
  
  RTVI <- 100 * (spectral850 - spectral725) - 10 * (spectral850 - spectral550_550)
  
  matVIs <- rbind(GRVI = GRVI,
                  NDVI = NDVI, 
                  GNDVI = GNDVI, 
                  BNDVI = BNDVI, 
                  NDRE = NDRE, 
                  CIgreen = CIgreen, 
                  GRVI = GRVI, 
                  CIRE = CIRE, 
                  MSR = MSR, 
                  MSRRE = MSRRE, 
                  NDI = NDI, 
                  BGVI = BGVI, 
                  EVI = EVI, 
                  VARI = VARI, 
                  ARI = ARI, 
                  mARI = mARI, 
                  BRVI = BRVI, 
                  BRENDI = BRENDI, 
                  GRENDI = GRENDI, 
                  RTVI = RTVI)
  matVIsSel <- matVIs[VIsName, ]
  plotMat <- rbind(plotMat0, matVIsSel)
  plotMat[is.infinite(plotMat)] <- 0
  return(plotMat)
}

SelectPlotCsvFromList <- function(allPlotCsvList, selectPlot, getOrRemove) {
  # selPlotCsvName <- paste0(selectPlot, ".csv")
  selPlotBind <- paste0("(", str_c(selectPlot, collapse = "|"), ")")
  selPlotPlace <- grep(selPlotBind, allPlotCsvList)
  if (getOrRemove == "get") {
    selPlotCsv <- allPlotCsvList[selPlotPlace]
  } else if (getOrRemove == "remove") {
    selPlotCsv <- allPlotCsvList[-selPlotPlace]
  }
  return(selPlotCsv)
}

####### 3.0.0.MTMFolder ###########
MakeCorMat <- function(cov_file, index_names) {
  #cov_file <- list.files(data_folder)[5]
  cov_df <- read.table(cov_file)
  cor_value <- cov2cor(as.matrix(cov_df))
  cor_value <- round(cor_value, digits = 2)
  colnames(cor_value) <- index_names
  rownames(cor_value) <- index_names
  return(cor_value)
}

####### 3.1.0.MTMprediction ###########
CsvToDf <- function(baseCsv, csvRowInd, csvColInd) {
  corMat <- as.matrix(baseCsv)
  colnames(corMat) <- csvColInd
  value <- c(corMat)
  colInd <- rep(csvColInd, each = nrow(corMat))
  rowInd <- rep(csvRowInd, ncol(corMat))
  longValue <- cbind(value, colInd, rowInd)
  longValue <- data.frame(longValue)
  longValue$value <- as.numeric(longValue$value)
  return(longValue)
}


