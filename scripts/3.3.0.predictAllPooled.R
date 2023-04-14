setwd("H:/TSMS_supple")
# machine <- "drone"
the_year <- "2019"
dataWeek <- c("week1", "week2", "week3", "week4", "week5", "week6")

library(RAINBOWR)
library(readr)
library(BGLR)
library(doParallel)
library(stringr)
library(ggplot2)
library(ggsci)
source("scripts/1.functionCode.R")
options(stringsAsFactors = FALSE)

# make the folder
predictFolder <- "2019_M100_Xacti4eyeCamera_Images/result/3.3.0.predictAllPooled"
if (!file.exists(predictFolder)) {
  dir.create(predictFolder)
}

# define the model
# "MS" is a kernel prediction model using MS kernel
# "G+MS" is a multi-kernel prediction model using G and MS kernel
# "G+Field" is a multi-kernel prediction model using G and Field kernel
# "NDVI" is a linear regression model using NDVI value
kernelName <- c("MS", "G+MS", "G+E", "G+E+GEI", "G+E+MS", "NDVI")
resultIndex <- c("correlation", "R2", "RMSE")

arrayResultRAINBOW <- array(NA, dim = c(length(dataWeek), 
                                        length(kernelName), length(resultIndex)))
dimnames(arrayResultRAINBOW) <- list(dataWeek, kernelName, resultIndex)


# save or read the seed
seedIndCsv <- "2019_M100_Xacti4eyeCamera_Images/result/seedInd.csv"
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(seedIndCsv, row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = seedIndCsv)
}

# read the genome data
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"

# data folder
phenoFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"

# set the field dist
fieldMap <- data.frame(read_csv("2019_M100_Xacti4eyeCamera_Images/field_data/2020_Tottori_Main_PlotPosition.csv"))
D2 <- as.matrix(dist(fieldMap[,c("x.real", "y.real")]))^2
# image(D2)
D2[1:100,101:200] <- Inf
D2[101:200,1:100] <- Inf
# image(D2)

D2.sub <- D2[1:100, 1:100]
h <- 1 / median(D2.sub[upper.tri(D2.sub)])
kPos <- exp(-h * D2)
# image(kPos)

# make field dist
dataRef <- read.csv(paste0(phenoFolder, "/week1_medianPheno.csv"), 
                    header = T, row.names = 1)
plotInd <- as.numeric(stringr::str_sub(rownames(dataRef), 4, 6))
kPosEach0 <- kPos[plotInd, plotInd]
rownames(kPosEach0) <- colnames(kPosEach0) <- rownames(dataRef)

cl <- makeCluster(10)
registerDoParallel(cl)
allResult <- foreach(dayInd = 1:length(dataWeek), .packages = c("BGLR", "RAINBOWR", "stringr")) %dopar% { 
  # dayInd <- 1
  the_day <- dataWeek[dayInd]
  
  saveFolder <- paste0(predictFolder, "/", the_day)
  if (!file.exists(saveFolder)) {
    dir.create(saveFolder)
  }
  dataAll0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"), 
                       header = T, row.names = 1)
  
  # remove "plot", and each spectral
  dataAll <- dataAll0[, -c(3, 7:11)]
  
  phenoDataAll0 <- na.omit(dataAll)
  phenoDataAll0 <- phenoDataAll0[phenoDataAll0$variety %in% colnames(amat0), ]
  
  # take log
  phenoDataAll0$dryWeight <- log(phenoDataAll0$dryWeight)
  phenoDataAll0$plantArea <- log(phenoDataAll0$plantArea)
  
  # meke treatment kernel
  treatmentMat <- model.matrix( ~ phenoDataAll0$treatment - 1)
  rownames(treatmentMat) <- rownames(phenoDataAll0)
  kernelTre <- tcrossprod(treatmentMat)
  # image(kernelTre)
  
  # dryweight
  dryWeight <- phenoDataAll0[, "dryWeight"]
  names(dryWeight) <- phenoDataAll0$variety
  
  #plant Area
  plantArea <- phenoDataAll0[, "plantArea"]
  names(plantArea) <- phenoDataAll0$variety
  
  # NDVI
  NDVI <- phenoDataAll0[, "NDVI"]
  names(NDVI) <- phenoDataAll0$variety
  
  # make MS data kernel
  phenoMs <- phenoDataAll0[, 6:ncol(phenoDataAll0)]
  phenoMs <- as.matrix(phenoMs)
  rownames(phenoMs) <- phenoDataAll0$variety
  phenoMsScaled <- scale(phenoMs)
  kernelMs <- tcrossprod(phenoMsScaled)
  
  # make amat
  amat <- amat0[names(dryWeight), names(dryWeight)]
  
  # make amat * treatment kernel
  kernelAmatTre <- amat * kernelTre
  # image(kernelAmatTre)
  
  # make amat list
  amatZ <- design.Z(pheno.labels = names(dryWeight),
                    geno.names = rownames(amat))
  
  amatList <- list(Z = amatZ,
                   K = amat)
  
  # make MS list
  msZ <- diag(1, ncol(kernelMs))
  rownames(msZ) <- colnames(msZ) <- rownames(kernelMs)
  
  msList <- list(Z = msZ,
                 K = kernelMs)
  
  # make G * treatment list
  amatTreZ <- design.Z(pheno.labels = names(dryWeight),
                       geno.names = rownames(kernelAmatTre))
  
  amatTreList <- list(Z = amatTreZ,
                      K = kernelAmatTre)
  
  # set each ZETA 
  ZETA_MS <- list(msList = msList)
  ZETA_GM <- list(amat = amatList,
                  msList = msList)
  ZETA_G <- list(amat = amatList)
  ZETA_GT <- list(amat = amatList,
                  amatTre = amatTreList)
  # make the fixed effect
  treatmentX0 <- treatmentMat
  
  NDVIX0 <- cbind(1, NDVI)
  colnames(NDVIX0) <- c("intercept", "NDVI")
  
  # make the fixed effect
  flowerDate <- phenoDataAll0$flower
  
  if (any(flowerDate == 0)) {
    flowerX0 <- cbind(1, flowerDate)
    rownames(flowerX0) <- names(dryWeight)
    multiX0 <- cbind(treatmentX0, flowerDate)
  } else {
    flowerX0 <- NULL
    multiX0 <- treatmentX0
  }
  rownames(multiX0) <- names(dryWeight)
  
  ZETAList <- list(list(ZETA = ZETA_MS,  X0 = flowerX0, nameInd = "MS"), 
                   list(ZETA = ZETA_GM,  X0 = flowerX0, nameInd = "G+MS"), 
                   list(ZETA = ZETA_G,   X0 = multiX0,  nameInd = "G+E"), 
                   list(ZETA = ZETA_GT,  X0 = multiX0,  nameInd = "G+E+GEI"), 
                   list(ZETA = ZETA_GM,  X0 = multiX0,  nameInd = "G+E+MS"), 
                   list(ZETA = ZETA_G,   X0 = NDVIX0,   nameInd = "NDVI"))
  
  
  eachResult <- lapply(ZETAList, function(eachList) {
    # eachList <- ZETAList[[4]]
    
    ZETA <- eachList$ZETA
    X0 <- eachList$X0
    nameInd <- eachList$nameInd
    
    # set the seed and cross-validation index
    rep5 <- rep(1:5, 1000)
    rep5 <- rep5[1:nrow(phenoDataAll0)]
    
    # make the case for the result
    resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
    colnames(resultEachSeed) <- resultIndex
    rownames(resultEachSeed) <- seedInd
    
    varietyName <- unique(phenoDataAll0[, "variety"])
    
    for (seedIndEach in 1:length(seedInd)) {
      # seedIndEach <- 1
      set.seed(seedInd[seedIndEach])
      varietyInd <- sample(rep5, length(varietyName), replace = F)
      # prepare the prediction box
      predictionDataRAINBOW <- rep(NA, nrow(phenoDataAll0))
      for (times in 1:5) {
        # times <- 1
        # remove the variety (which tested) from training data
        testVariety <- varietyName[varietyInd == times]
        testInd <- phenoDataAll0$variety %in% testVariety
        
        dryWeightNa <- dryWeight
        dryWeightNa[testInd] <- NA
        
        if (ZETA == "lm") {
          
          lmFit <- lm(dryWeight[!testInd] ~ X0[!testInd, "NDVI"])
          a <- lmFit$coefficients[2]
          b <- lmFit$coefficients[1]
          lmPred <- a * X0[testInd, "NDVI"] + b
          
          predictionDataRAINBOW[testInd] <- lmPred
        } else {
          # use RAINBOWR
          resEM3 <- EM3.cpp(y = dryWeightNa,
                            X0 = X0, n.core = 1,
                            ZETA = ZETA)
          
          # input the predicted data
          predictionDataRAINBOW[testInd] <- resEM3$y.pred[testInd]
        }
      }
      predictData <- (predictionDataRAINBOW)
      obsData <- (dryWeight)
      # calculate the R2 and RMSE
      correlation <- cor(obsData, predictData)
      R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
      RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
      
      # input the result
      resultEachSeed[seedIndEach, ] <- c(correlation, R2, RMSE)
      
      # make the plot
      xlim <- ylim <- range(predictData, obsData)
      
      png(paste0(saveFolder, "/predictionPlot",
                 seedIndEach, "_", nameInd, ".png"))
      plot((obsData), (predictData),
           xlim = xlim, ylim = ylim,
           main = paste0("RAINBOW prediction ", nameInd, " r = ", round(correlation, 2)))
      abline(0, 1, col = 2, lty = 2)
      dev.off()
    }
    return(apply(resultEachSeed, 2, mean))
  })
  eachDayResult <- do.call(what = rbind, args = eachResult)
  rownames(eachDayResult) <- kernelName
  write.csv(eachDayResult, file = paste0(saveFolder, "/eachDayResult.csv"))
  return(eachDayResult)
}
stopCluster(cl)



########## visualize the result #############
resultFolder <- predictFolder

for (eachDay in dataWeek) {
  # eachDay <- dataWeek[2]
  resultFolderEach <- paste0(resultFolder, "/", eachDay)
  resultFile <- list.files(resultFolderEach, pattern = ".csv", full.names = T)
  resultDf <- read.csv(resultFile[[1]], header = T, row.names = 1)
  arrayResultRAINBOW[eachDay, , ] <- as.matrix(resultDf)
}
arrayResultRAINBOW <- arrayResultRAINBOW[, -3, ]

ShowArrayResult <- function(resultArray, resultArrayInd, resultName) {
  eachResultMat <- resultArray[, , resultArrayInd]
  eachResultDf <- CsvToDf(baseCsv = eachResultMat, 
                          csvRowInd = rownames(eachResultMat), 
                          csvColInd = colnames(eachResultMat))
  colnames(eachResultDf) <- c("value", "model", "Week")
  eachResultDf$Week <- as.factor(eachResultDf$Week)
  eachResultDf$Week <- as.numeric(eachResultDf$Week)
  
  g <- ggplot(eachResultDf, aes(x = Week, y = value, colour = model)) + 
    geom_line(size = 1) + 
    scale_color_manual(values = c("steelblue2", "coral4", "springgreen3", "tan1", "gray40")) + 
    scale_x_continuous(breaks = c(1:6), labels = c(1:6))
  return(g)
}

gCorrelation <- ShowArrayResult(resultArray = arrayResultRAINBOW, resultArrayInd = 1, resultName = "correlation")
gR2 <- ShowArrayResult(resultArray = arrayResultRAINBOW, resultArrayInd = 2, resultName = "R2")
gRMSE <- ShowArrayResult(resultArray = arrayResultRAINBOW, resultArrayInd = 3, resultName = "RMSE")

# write down as .pdf
pdf(paste0(predictFolder, "/resultLinePooled.pdf"))
print(gCorrelation)
print(gR2)
print(gRMSE)
dev.off()

# write down as .png
png(paste0(predictFolder, "/resultLineCorrelationPooled.png"), 
    height = 1440, width = 1440, res = 216)
print(gCorrelation)
dev.off()

png(paste0(predictFolder, "/resultLineR2Pooled.png"), 
    height = 1440, width = 1440, res = 216)
print(gR2)
dev.off()

png(paste0(predictFolder, "/resultLineRMSEPooled.png"), 
    height = 1440, width = 1440, res = 216)
print(gRMSE)
dev.off()