
# machine <- "drone"
the_year <- "2019"
dataWeek <- c("week1", "week2", "week3", "week4", "week5", "week6")

#' # 0. Read packages
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
kernelFolder <- "2019_M100_Xacti4eyeCamera_Images/result/3.1.2.kernelPredictionWithFlowerDate"
if (!dir.exists(kernelFolder)) {
  dir.create(kernelFolder)
}

# data folder
phenoFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"

# set the condition
condition <- c("W1", "W2", "W3", "W4")

# define the model "G" is a genomic prediction model
# "MS" is a kernel prediction model using MS kernel
# "G+MS" is a multi-kernel prediction model using G and MS kernel
# "G+Field" is a multi-kernel prediction model using G and Field kernel
# "NDVI" is a linear regression model using NDVI value
kernelName <- c("G", "MS", "G+MS", "G+Field", "NDVI")
resultIndex <- c("Correlation", "R2", "RMSE")

# make the array to input correlation
arrayResultRAINBOW <- array(NA, dim = c(length(condition), length(kernelName), length(resultIndex)))
dimnames(arrayResultRAINBOW) <- list(condition, kernelName, resultIndex)


# save or read the seed
seedIndCsv <- "2019_M100_Xacti4eyeCamera_Images/result/seedInd.csv"
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(seedIndCsv, row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = seedIndCsv)
}

# read the field ID and plot
plotIdTable <- read.csv("2019_M100_Xacti4eyeCamera_Images/field_data/PlotIDtable_2019_Jul.csv", 
                        row.names = 1)
plotIdGeno <- plotIdTable$Name

# read the genome data
amat0 <- as.matrix(read.csv("genome/amat173583SNP.csv", row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"

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

cl <- makeCluster(7)
registerDoParallel(cl)
allResult <- foreach(dayInd = 1:length(dataWeek), .packages = c("BGLR", "RAINBOWR", "stringr")) %dopar% { 
  # dayInd <- 1
  the_day <- dataWeek[dayInd]
  
  saveFolder <- paste0(kernelFolder, "/", the_day)
  if (!file.exists(saveFolder)) {
    dir.create(saveFolder)
  }
  # phenotype data
  phenoEachDay0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"),
                            header = TRUE, row.names = 1)
  dim(phenoEachDay0)
  # remove "plot",  "flower", and each spectral
  phenoEachDay <- phenoEachDay0[, -c(3, 7:11)]
  
  phenoDataAll0 <- na.omit(phenoEachDay)
  phenoDataAll0 <- phenoDataAll0[phenoDataAll0$variety %in% colnames(amat0), ]
  # plot(phenoDataAll0$dryWeight, phenoDataAll0$plantArea)
  # 
  phenoDataAll0$dryWeight <- log(phenoDataAll0$dryWeight)
  phenoDataAll0$plantArea <- log(phenoDataAll0$plantArea)
  for (conditionEach in condition) {
    # conditionEach <- condition[1]
    
    # select the training and test condition data
    phenoDataAll <- phenoDataAll0[phenoDataAll0$treatment == conditionEach, ]
    dim(phenoDataAll)

    # dryweight
    dryWeight0 <- phenoDataAll[, "dryWeight"]
    dryWeight <- scale(phenoDataAll[, "dryWeight"])
    names(dryWeight) <- phenoDataAll$variety
    # image(tcrossprod(scale(dryWeight)))
    
    # NDVI
    NDVI0 <- phenoDataAll[, "NDVI"]
    NDVI <- scale(phenoDataAll[, "NDVI"])
    names(NDVI) <- phenoDataAll$variety
    
    # make MS data kernel
    phenoMs <- phenoDataAll[, 6:ncol(phenoDataAll)]
    phenoMs <- as.matrix(phenoMs)
    rownames(phenoMs) <- phenoDataAll$variety
    phenoMsScaled <- scale(phenoMs)
    kernelMs <- tcrossprod(phenoMsScaled)
    
    # make amat
    amat <- amat0[names(dryWeight), names(dryWeight)]

    # make field dist
    kPosEach0 <- kPos
    plotInd <- as.numeric(stringr::str_sub(rownames(phenoDataAll), 4, 6))
    kPosEach <- kPosEach0[plotInd, plotInd]
    # dim(kPosEach)
    rownames(kPosEach) <- colnames(kPosEach) <- phenoDataAll$variety
    
    
    # make amat list
    amatZ <- design.Z(pheno.labels = names(dryWeight),
                      geno.names = rownames(amat))
    
    amatList <- list(Z = amatZ,
                     K = amat)
    
    # make MS list
    msZ <- design.Z(pheno.labels = names(dryWeight),
                    geno.names = rownames(kernelMs))
    # msZ <- diag(1, ncol(kernelMs))
    
    msList <- list(Z = msZ,
                   K = kernelMs)
    
    # make filed kernel list
    fieldZ <- design.Z(pheno.labels = names(dryWeight),
                       geno.names = rownames(kPosEach))
    fieldList <- list(Z = fieldZ,
                      K = kPosEach)
    
    # set each ZETA 
    ZETA_G <- list(amat = amatList)
    ZETA_MS <- list(msList = msList)
    ZETA_GM <- list(amat = amatList,
                    msList = msList)
    ZETA_GF <- list(amat = amatList, 
                    fieldList = fieldList)
    
    # make the fixed effect
    flowerDate <- phenoDataAll$flower
    if (any(flowerDate == 0)) {
      flowerX0 <- cbind(1, flowerDate)
      rownames(flowerX0) <- names(dryWeight)
      colnames(flowerX0) <- c("intercept", "flowering")
    } else {
      flowerX0 <- NULL
    }
    
    # set the fixed effect
    NDVIX0 <- cbind(1, NDVI)
    colnames(NDVIX0) <- c("intercept", "NDVI")
    
    
    ZETAList <- list(list(ZETA = ZETA_G,  X0 = flowerX0, nameInd = "G"), 
                     list(ZETA = ZETA_MS, X0 = flowerX0, nameInd = "MS"), 
                     list(ZETA = ZETA_GM, X0 = flowerX0, nameInd = "G+MS"), 
                     list(ZETA = ZETA_GF, X0 = flowerX0, nameInd = "G+Field"), 
                     list(ZETA = "lm", X0 = NDVIX0, nameInd = "NDVI"))
    
    eachResult <- lapply(ZETAList, function(eachList) {
      # eachList <- ZETAList[[5]]
      
      ZETA <- eachList$ZETA
      X0 <- eachList$X0
      nameInd <- eachList$nameInd
      
      # set the seed and cross-validation index
      rep5 <- rep(1:5, 1000)
      rep5 <- rep5[1:nrow(phenoDataAll)]

      # make the case for the result
      resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
      colnames(resultEachSeed) <- resultIndex
      rownames(resultEachSeed) <- seedInd
      
      
      for (seedIndEach in 1:length(seedInd)) {
        # seedIndEach <- 1
        set.seed(seedInd[seedIndEach])
        crossCvInd <- sample(rep5, nrow(phenoDataAll), replace = F)
        
        # prepare the prediction box
        predictionDataRAINBOW <- rep(NA, nrow(phenoDataAll))
        for (times in 1:5) {
          # times <- 1
          testInd <- crossCvInd == times
          
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
        # resEM3$weights
        predictData <- (predictionDataRAINBOW)
        obsData <- (dryWeight)
        # calculate the R2 and RMSE
        correlation <- cor(obsData, predictData)
        R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
        RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
        
        # input the result
        resultEachSeed[seedIndEach, ] <- c(correlation, R2, RMSE)
        # arrayResultEachDay[ZETA, ] <- c(correlation, R2, RMSE)
        
        # make the plot
        xlim <- ylim <- range(predictData, obsData)
        
        png(paste0(saveFolder, "/", conditionEach, "_predictionPlot",
                   seedIndEach, "_", nameInd, ".png"))
        plot((obsData), (predictData), 
             xlim = xlim, ylim = ylim, 
             main = paste0("RAINBOW prediction ", nameInd, " r = ", round(correlation, 2)))
        abline(0, 1, col = 2, lty = 2)
        dev.off()
      }
      return(apply(resultEachSeed, 2, mean))
    })
    
    eachConditionResult <- do.call(what = rbind, args = eachResult)
    rownames(eachConditionResult) <- kernelName
    write.csv(eachConditionResult, file = paste0(saveFolder, "/", conditionEach, "_Result.csv"))
  }
}
stopCluster(cl)


allDayResultDfList <- lapply(dataWeek, function(eachDay) {
  # eachDay <- dataWeek[1]
  resultFolder <- paste0(kernelFolder, "/", eachDay)
  eachDayResultList <- list.files(resultFolder, pattern = ".csv", full.names = T)
  eachDayResultDfList <- lapply(eachDayResultList, function(eachDayEachCondition) {
    # eachDayEachCondition <- eachDayResultList[[1]]
    eachCondition0 <- str_split(eachDayEachCondition, pattern = "/")[[1]][5]
    eachCondition <- str_sub(eachCondition0, 1, 2)
    eachDayEachConditionCsv <- read.csv(eachDayEachCondition, header = T, row.names = 1)
    eachDayEachConditionDf0 <- CsvToDf(baseCsv = eachDayEachConditionCsv, 
                                       csvRowInd = rownames(eachDayEachConditionCsv), 
                                       csvColInd = colnames(eachDayEachConditionCsv))
    eachDayEachConditionDf <- cbind(eachDayEachConditionDf0, condition = eachCondition)
    return(eachDayEachConditionDf)
  })
  eachDayResultDf0 <- do.call(what = rbind, args = eachDayResultDfList)
  eachDayResultDf <- cbind(eachDayResultDf0, week = eachDay)
  return(eachDayResultDf)
})

allDayResultDf0 <- do.call(what = rbind, args = allDayResultDfList)
write.csv(allDayResultDf0, paste0(kernelFolder, "/resultLongDf.csv"))
colnames(allDayResultDf0) <- c("value", "index", "model", "condition", "week")
allDayResultDf <- allDayResultDf0

allDayResultDf$week <- as.factor(allDayResultDf$week)
allDayResultDf$week <- as.numeric(allDayResultDf$week)
resultIndList <- list("Correlation", "R2", "RMSE")

# set the dataframe for ggplot
resultForGG <- allDayResultDf
colnames(resultForGG)[c(3, 5)] <- c("Model", "Week")
conditionList <- list(list("W1", "C"), 
                      list("W2", "W5"), 
                      list("W3", "W10"), 
                      list("W4", "D"))
for (eachCondition in conditionList) {
  # eachCondition <- conditionList[[1]]
  conditionOld <- eachCondition[[1]]
  conditionNew <- eachCondition[[2]]
  resultForGG$condition[resultForGG$condition == conditionOld] <- conditionNew
}
resultForGG$condition <- factor(resultForGG$condition, 
                                levels = c("C", "W5", "W10", "D"))

ggLineList <- lapply(resultIndList, function(eachResultInd) {
  # eachResultInd <- resultIndList[[1]]
  allDayResultEach <- resultForGG[resultForGG$index == eachResultInd, ]
  xmin <- min(allDayResultEach$Week)
  xmax <- max(allDayResultEach$Week)
  g <- ggplot(allDayResultEach, aes(x = Week, y = value, colour = Model, linetype = Model)) + 
    facet_wrap(~ condition) + 
    geom_line(size = 1) + 
    # ylim(ymin, ymax) +
    scale_linetype_manual(values = c("solid", "dashed", "solid", "solid", "solid")) + 
    scale_color_manual(values = c("steelblue2", "ivory4", "springgreen3", "tan1", "gray40")) + 
    labs(y = paste0(eachResultInd)) + 
    scale_x_continuous(breaks = seq(xmin, xmax, by = 1),
                       labels = seq(xmin, xmax, by = 1))
  return(g)
})

pdf(paste0(kernelFolder, "/resultLine.pdf"))
print(ggLineList[[1]])
print(ggLineList[[2]])
print(ggLineList[[3]])
dev.off()

# write down as .png
png(paste0(kernelFolder, "/resultLineCorrelation.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineList[[1]])
dev.off()

png(paste0(kernelFolder, "/resultLineR2.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineList[[2]])
dev.off()

png(paste0(kernelFolder, "/resultLineRMSE.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineList[[3]])
dev.off()

# remove the kernel of "G+Field"
allDayResultDfSel <- resultForGG[resultForGG$Model != "G+Field", ]
ggLineListSel <- lapply(resultIndList, function(eachResultInd) {
  # eachResultInd <- resultIndList[[1]]
  allDayResultEach <- allDayResultDfSel[allDayResultDfSel$index == eachResultInd, ]
  xmin <- min(allDayResultEach$Week)
  xmax <- max(allDayResultEach$Week)
  g <- ggplot(allDayResultEach, aes(x = Week, y = value, colour = Model, linetype = Model)) + 
    facet_wrap(~ condition) + 
    geom_line(size = 1) + 
    # ylim(ymin, ymax) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "dashed")) + 
    scale_color_manual(values = c("steelblue2", "springgreen3", "tan1", "tomato2")) + 
    labs(y = paste0(eachResultInd)) + 
    scale_x_continuous(breaks = seq(xmin, xmax, by = 1),
                       labels = seq(xmin, xmax, by = 1))
  return(g)
})

# write down as .png
png(paste0(kernelFolder, "/resultLineCorrelationSel.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineListSel[[1]])
dev.off()

png(paste0(kernelFolder, "/resultLineR2Sel.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineListSel[[2]])
dev.off()

png(paste0(kernelFolder, "/resultLineRMSESel.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineListSel[[3]])
dev.off()
