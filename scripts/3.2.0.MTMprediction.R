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

# make the folder
mtmFolder <- "2019_M100_Xacti4eyeCamera_Images/result/3.2.0.MTMprediction"
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

# set the condition
condition <- c("W1", "W2", "W3", "W4")

# make the data folder
phenoFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"

# read the amat
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
# amat <- as.matrix(amat0)
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"


#### calculate the MTM ######
cl <- makeCluster(10)
registerDoParallel(cl)
resultAll <- foreach(dayInd = 1:length(dataWeek), .packages = c("MTM", "corrplot")) %dopar% {
  # the_day <- dataWeek[[4]]
  the_day <- dataWeek[[dayInd]]
  
  # read the phenotypic data and extract the data
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
  
  # make the matrix to input the result
  corEachDaySeed <- matrix(NA, nrow = length(condition), ncol = length(seedInd))
  rownames(corEachDaySeed) <- condition
  colnames(corEachDaySeed) <- seedInd
  
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
    
    # extract the data of the target treatment
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
    
    dryWeight <- MTM_Y[, "dryWeight"]
    
    # make the vector to input the correlation between observed value and predicted value
    corMTM <- rep(NA, length(seedInd))
    
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
          res_MTM <- MTM_2(Y = MTM_Y_Na, K = MTM_K, resCov = MTM_resCov,
                           nIter = 12000, burnIn = 2000, thin = 20,
                           saveAt = paste0(new_folder, "/", treatment_i, '_', seedIndEach, "_", times, "_"))
        )
        
        # correlation between prediction & true value
        predictionData[crossCvInd == times] <- res_MTM$YHat[, "dryWeight"][crossCvInd == times]
      }
      corMTM[seedIndEach] <- cor(predictionData, dryWeight)
      
      xlim <- ylim <- range(predictionData, dryWeight)
      png(paste0(new_folder, "/predictionPlot", 
                 treatment_i, "_", seedIndEach, ".png"))
      plot(dryWeight, predictionData, 
           main = paste0(" MTM prediction ", treatment_i), 
           xlim = xlim, ylim = ylim)
      abline(0, 1, col = 2, lty = 2)
      dev.off()
    }
    write.csv(corMTM, paste0(new_folder, "/", treatment_i, "_result.csv"))
    corEachDaySeed[treatment_i, ] <- corMTM
  }
}

stopCluster(cl)

######## read the each day data ############
allDayResultDfList <- lapply(dataWeek, function(eachDay) {
  # eachDay <- dataWeek[1]
  resultFolder <- paste0(mtmFolder, "/", the_year, eachDay)
  eachDayResultList <- list.files(resultFolder, pattern = "_result.csv", full.names = T)
  eachDayResultDfList <- lapply(eachDayResultList, function(eachDayEachCondition) {
    # eachDayEachCondition <- eachDayResultList[[1]]
    eachCondition0 <- str_split(eachDayEachCondition, pattern = "/")[[1]][5]
    eachCondition <- str_sub(eachCondition0, 1, 2)
    eachDayEachConditionCsv <- read.csv(eachDayEachCondition, header = T, row.names = 1)
    return(mean(as.matrix(eachDayEachConditionCsv)))
  })
  eachDayResultDf0 <- do.call(what = rbind, args = eachDayResultDfList)
  return(eachDayResultDf0)
})

allDayResultDf <- do.call(what = cbind, args = allDayResultDfList)
rownames(allDayResultDf) <- condition
colnames(allDayResultDf) <- dataWeek
write.csv(allDayResultDf, paste0(mtmFolder, "/result.csv"))


########## visualize the result #########
resultCsvMT <- paste0(mtmFolder, "/result.csv")
resultMT <- read.csv(resultCsvMT, header = T, row.names = 1)
resultMTlong0 <- CsvToDf(baseCsv = resultMT, 
                         csvRowInd = rownames(resultMT), 
                         csvColInd = colnames(resultMT))
resultMTlong <- cbind(resultMTlong0, "MT")
colnames(resultMTlong) <- c("value", "week", "treatment", "model")
write.csv(resultMTlong, paste0(mtmFolder, "/resultLongDf.csv"))

resultCsvST <- "2019_M100_Xacti4eyeCamera_Images/result/3.1.1.kernelPrediction/resultLongDf.csv"
resultST0 <- read.csv(resultCsvST, header = T, row.names = 1)
resultSTlong0 <- resultST0[resultST0$colInd == "Correlation" & resultST0$rowInd == "G", ]
resultSTlong <- data.frame(resultSTlong0$value, resultSTlong0$week, resultSTlong0$condition, "G")
colnames(resultSTlong) <- c("value", "week", "treatment", "model")

resultAllLong <- rbind(resultMTlong, resultSTlong)
resultAllLong$week <- as.factor(resultAllLong$week)
resultAllLong$week <- as.numeric(resultAllLong$week)
g <- ggplot(resultAllLong, aes(x = week, y = value, colour = model)) + 
  facet_wrap(~ treatment) + 
  geom_line(size = 1) + 
  # ylim(ymin, ymax) +
  scale_color_manual(values = c("steelblue2", "tomato2")) + 
  labs(title = "MT and ST model prediction",  
       y = "correlation")

# save the figure
png(paste0(mtmFolder, "/predictResultLine.png"), 
    height = 1440, width = 1440, res = 216)
print(g)
dev.off()

