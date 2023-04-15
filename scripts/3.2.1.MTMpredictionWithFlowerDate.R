
# machine <- "drone"
the_year <- "2019"
dataWeek <- c("week1", "week2", "week3", "week4", "week5", "week6")

#' # 0. Read packages
library(rrBLUP)
library(MASS)
library(tidyr)
library(MTM)
library(corrplot)
library(outliers)
library(stringr)
library(reshape2)
library(date)
library(ggplot2)
library(ggsci)
library(doParallel)
source("scripts/1.functionCode.R")
source("scripts/1.MTM_2.R")
options(stringsAsFactors = FALSE)

# make a folder
mtmFolder <- "2019_M100_Xacti4eyeCamera_Images/result/3.2.1.MTMpredictionWithFlowerDate"
if (!file.exists(mtmFolder)) {
  dir.create(mtmFolder)
}

# save or read the seed
seedIndCsv <- "2019_M100_Xacti4eyeCamera_Images/result/seedInd.csv"
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(seedIndCsv, row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = seedIndCsv)
}

# data folder
phenoFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"

# set the condition
condition <- c("W1", "W2", "W3", "W4")
resultIndex <- c("Correlation", "R2", "RMSE")

# make the array to input correlation
arrayResultMTM <- array(NA, dim = c(length(condition), length(resultIndex)))
dimnames(arrayResultMTM) <- list(condition, resultIndex)

# read the amat
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"

#### calculate the MTM ######
if (length(list.files(mtmFolder)) < 7) {
  cl <- makeCluster(10)
  registerDoParallel(cl)
  resultAll <- foreach(dayInd = 1:length(dataWeek), .packages = c("MTM", "corrplot")) %dopar% {
    # the_day <- dataWeek[[1]]
    the_day <- dataWeek[[dayInd]]
    
    # read the phenotypic data
    dataFrame0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"),
                           header = TRUE, row.names = 1)
    df <- dataFrame0[, c("variety", "treatment", "flower", "plantArea", "dryWeight", 
                         "GRVI",
                         "NDVI", 
                         "NDRE", 
                         "NDI", 
                         "RTVI")]
    df <- na.omit(df)
    df <- df[, c("variety", "treatment", "flower", "dryWeight", 
                 "GRVI",
                 "NDVI", 
                 "NDRE", 
                 "NDI", 
                 "RTVI")]
    
    # set the treatments and traits
    n_treatment <- length(unique(df$treatment))
    name_treatment <- c("W1", "W2", "W3", "W4")
    name_trait <- colnames(df[, c(4:ncol(df))])
    n_trait <- length(name_trait)
    
    phenodat <- df
    phenodat$dryWeight <- log(df$dryWeight)
    colnames(phenodat)[c(1, 2)] <- c("line", "treatment")
    # head(phenodat)
    
    
    # 3. Estimate genotypic values with MTM
    # Estimate genotypic values for each condition,
    new_folder <- paste0(mtmFolder, "/", the_year, the_day)
    if (!file.exists(new_folder)) {
      dir.create(new_folder)
    }
    
    for (treatment_i in name_treatment) {
      # treatment_i <- name_treatment[1]
      
      # choose data of the target treatment
      phenodat_i <- phenodat[phenodat$treatment == treatment_i, ]
      phenodat_i <- na.omit(phenodat_i)
      
      phenodat_i <- phenodat_i[phenodat_i$line %in% colnames(amat0), ]
      amat <- amat0[as.vector(phenodat_i$line), as.vector(phenodat_i$line)]
      name_line <- rownames(amat)
      
      # Y
      MTM_Y <- as.matrix(phenodat_i[, -c(1:3)])
      MTM_Y <- scale(MTM_Y)
      rownames(MTM_Y) <- phenodat_i$line
      
      # structure of the variance-covariance matrix of genotypic values
      MTM_K <- list(list(
        K = amat,
        COV = list(type = 'UN',
                   df0 = n_trait,
                   S0 = diag(n_trait))
      ))
      
      # structure of the variance-covariance matrix of residual effects
      MTM_resCov <- list(type = 'UN',
                         df0 = n_trait,
                         S0 = diag(n_trait))
      
      # structure of fixed effect
      flower <- phenodat_i$flower
      if (all(flower == 1)) {
        MTM_XF <- NULL
      } else {
        flowerFactor <- as.factor(flower)
        MTM_XF <- model.matrix( ~ flowerFactor - 1)
        MTM_XF <- matrix(MTM_XF[, 1], ncol = 1)
      }
      
      
      dryWeight <- MTM_Y[, "dryWeight"]
      
      # make the array for the result
      resultEachSeed <- matrix(NA, ncol = length(resultIndex), nrow = length(seedInd))
      colnames(resultEachSeed) <- resultIndex
      rownames(resultEachSeed) <- seedInd
      
      rep5 <- rep(1:5, length(dryWeight))[1:length(dryWeight)]
      
      for (seedIndEach in 1:length(seedInd)) {
        # seedIndEach <- 1
        set.seed(seedInd[seedIndEach])
        crossCvInd <- sample(rep5, length(dryWeight), replace = F)
        predictionData <- rep(NA, length(dryWeight))
        
        for (times in 1:5) {
          # times <- 1
          
          # enter NA
          MTM_Y_Na <- MTM_Y
          MTM_Y_Na[crossCvInd == times, "dryWeight"] <- NA
          
          # estimation
          MTM_output <- capture.output(
            res_MTM <- MTM_2(XF = MTM_XF, Y = MTM_Y_Na, K = MTM_K, resCov = MTM_resCov,
                             nIter = 48000, burnIn = 8000, thin = 40,
                             saveAt = paste0(new_folder, "/", treatment_i, '_', seedIndEach, "_", times, "_"))
          )

          # put into the predicted value
          if (is.null(MTM_XF)) {
            predictEach <- (res_MTM$mu + res_MTM$K[[1]]$U)[, 1]
          } else {
            predictEach <- (res_MTM$mu + MTM_XF %*% res_MTM$B.f + res_MTM$K[[1]]$U)[, 1]
          }
          predictionData[crossCvInd == times] <- predictEach[crossCvInd == times]
        }
        
        predictData <- predictionData
        obsData <- dryWeight
        
        # calculate the R2 and RMSE
        correlation <- cor(obsData, predictData)
        R2 <- 1 - sum((obsData - predictData) ^ 2) / sum((obsData - mean(obsData)) ^ 2)
        RMSE <- sqrt(sum((obsData - predictData) ^ 2) / length(obsData))
        
        # input the result
        resultEachSeed[seedIndEach, ] <- c(correlation, R2, RMSE)
        
        xlim <- ylim <- range(predictData, obsData)
        png(paste0(new_folder, "/predictionPlot", 
                   treatment_i, "_", seedIndEach, ".png"))
        plot(obsData, predictData, 
             main = paste0(" MTM prediction ", treatment_i), 
             xlim = xlim, ylim = ylim)
        abline(0, 1, col = 2, lty = 2)
        dev.off()
      }
      write.csv(resultEachSeed, paste0(new_folder, "/", treatment_i, "_AllResult.csv"))
      resultEachCondition <- as.matrix(apply(resultEachSeed, 2, mean))
      resultEachCondition <- t(resultEachCondition)
      resultEachCondition <- cbind(resultEachCondition, sd = sd(resultEachSeed[, 1]))
      rownames(resultEachCondition) <- "MT"
      write.csv(resultEachCondition, paste0(new_folder, "/", treatment_i, "_result.csv"))
    }
  }
  
  stopCluster(cl)
}

allDayResultDfList <- lapply(dataWeek, function(eachDay) {
  # eachDay <- dataWeek[1]
  resultFolder <- paste0(mtmFolder, "/", the_year, eachDay)
  eachDayResultList <- list.files(resultFolder, pattern = "_result.csv", full.names = T)
  eachDayResultDfList <- lapply(eachDayResultList, function(eachDayEachCondition) {
    # eachDayEachCondition <- eachDayResultList[[1]]
    # eachCondition0 <- str_split(eachDayEachCondition, pattern = "/")[[1]][7]
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
write.csv(allDayResultDf0, paste0(mtmFolder, "/resultLongDf.csv"))
colnames(allDayResultDf0) <- c("value", "index", "model", "condition", "week")
allDayResultDf <- allDayResultDf0



########## visualize the result #########
resultCsvMT <- paste0(mtmFolder, "/resultLongDf.csv")
resultMTlong <- read.csv(resultCsvMT, header = T, row.names = 1)

resultCsvST <- "2019_M100_Xacti4eyeCamera_Images/result/3.1.2.kernelPredictionWithFlowerDate/resultLongDf.csv"
resultST0 <- read.csv(resultCsvST, header = T, row.names = 1)
resultSTlong <- resultST0[resultST0$rowInd == "G", ]

resultAllLong <- rbind(resultMTlong, resultSTlong)
colnames(resultAllLong) <- c("value", "index", "model", "condition", "week")

resultAllLong$week <- as.factor(resultAllLong$week)
resultAllLong$week <- as.numeric(resultAllLong$week)

resultIndList <- list("Correlation", "R2", "RMSE")

# set the dataframe for ggplot
resultForGG <- resultAllLong
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
  xmin <- min(resultForGG$Week)
  xmax <- max(resultForGG$Week)
  
  g <- ggplot(allDayResultEach, aes(x = Week, y = value, colour = Model, linetype = Model)) +
    facet_wrap(~ condition) +
    geom_line(size = 1) +
    # ylim(ymin, ymax) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
    scale_color_manual(values = c("steelblue2", "tomato2")) +
    labs(y = paste0(eachResultInd)) + 
    scale_x_continuous(breaks = seq(xmin, xmax, by = 1),
                       labels = seq(xmin, xmax, by = 1))
  return(g)
})

# write down as .png
png(paste0(mtmFolder, "/resultLineCorrelation.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineList[[1]])
dev.off()

png(paste0(mtmFolder, "/resultLineR2.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineList[[2]])
dev.off()

png(paste0(mtmFolder, "/resultLineRMSE.png"), 
    height = 1440, width = 1440, res = 216)
print(ggLineList[[3]])
dev.off()
