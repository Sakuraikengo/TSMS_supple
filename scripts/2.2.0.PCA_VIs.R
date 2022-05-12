#####prcomp######
# do PCA for the each day
machine <- "drone"
the_year <- "2019"
dataWeek <- c("week1", "week2", "week3", "week4", "week5", "week6")

#' # 0. Read packages
library(BGLR)
library(doParallel)
# library(gaston)
source("R/1.functionCode.R")
options(stringsAsFactors = FALSE)
library(lme4)
library(RAINBOWR)
library(stringr)
library(psych)
library(ggplot2)
library(corrplot)

# make the save folder
PCA_folder <- "2019_M100_Xacti4eyeCamera_Images/result/2.2.0.PCA_VIs"
if (!file.exists(PCA_folder)) {
  dir.create(PCA_folder)
}

# read the genome data
amat0 <- as.matrix(read.csv("genome/amat173583SNP.csv", row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"


# data folder
phenoFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"

# matrix to input the PCA result (Cumulative proportion)
matCP <- matrix(NA, nrow = length(dataWeek), ncol = 5)
rownames(matCP) <- dataWeek
colnames(matCP) <- paste0("PC", 1:5)


# array to input the variance of MS data
VIsName <- c("GRVI", "NDVI", "NDRE", "NDI", "RTVI")
condition <- c("W1", "W2", "W3", "W4")

arrayCor <- array(NA, dim = c(length(dataWeek), 
                              length(VIsName), length(condition)))
dimnames(arrayCor) <- list(dataWeek, VIsName, condition)


for (dayInd in 1:length(dataWeek)) {
  # dayInd <- 1
  the_day <- dataWeek[dayInd]
  dataAll0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"), 
                       header = T, row.names = 1)
  dataAll <- dataAll0[, -c(3, 4, 7:11)]
  dataAll <- na.omit(dataAll)
  dataAll <- dataAll[dataAll$variety %in% colnames(amat0), ]
  
  propMat <- matrix(NA, nrow = 4, ncol = 5)
  for (eachCondition in unique(dataAll$treatment)) {
    # eachCondition <- unique(dataAll$treatment)[1]
    dataEach <- dataAll[dataAll$treatment == eachCondition, ]
    dim(dataEach)
    msEach <- dataEach[, 5:ncol(dataEach)]
    
    resEach <- prcomp(msEach, scale = T)
    resSum <- summary(resEach)
    
    # correlation among AGB and PCs
    vecCor <- cor(cbind(dataEach$dryWeight, resEach$x))[1, 2:6]
    if (vecCor[1] < 0) {
      vecCor <- -(vecCor)
    }
    arrayCor[dayInd, , eachCondition] <- vecCor
    
    rownames(propMat) <- unique(dataAll$treatment)
    colnames(propMat) <- colnames(resSum$importance)
    propMat[eachCondition, ] <- resSum$importance[3, ]
    
    png(paste0(PCA_folder, "/pcaBarPlot_", the_day, "_", eachCondition, ".png"), 
        height = 1440, width = 1440, res = 216)
    
    plot(resEach)
    
    dev.off()
    
    png(paste0(PCA_folder, "/pcaBiPlot_", the_day, "_", eachCondition, ".png"), 
        height = 1440, width = 1440, res = 216)
    op <- par(mfrow = c(2,2))
    biplot(resEach, choices = 1:2)
    biplot(resEach, choices = 3:4)
    
    factor.loadings <- cor(msEach, resEach$x[,1:4])
    #factor.loadings
    
    theta <- 2 * pi * (0:100 / 100)
    x <- cos(theta)
    y <- sin(theta)
    plot(factor.loadings[,1:2], xlim = c(-1,1), ylim = c(-1,1), pch = " ")
    text(factor.loadings[,1:2], rownames(factor.loadings), col = "red")
    lines(x, y, col = "gray")
    abline(v = 0, h = 0)
    plot(factor.loadings[,3:4], xlim = c(-1,1), ylim = c(-1,1), pch = " ")
    text(factor.loadings[,3:4], rownames(factor.loadings), col = "red")
    lines(x, y, col = "gray")
    abline(v = 0, h = 0)
    par(op)
    
    dev.off()
  }
  write.csv(round(propMat, 2), 
            paste0(PCA_folder, "/", the_day, "_importance.csv"))
  
  
  phenoData <- dataAll[, 3:ncol(dataAll)]
  msData <- dataAll[, 5:ncol(dataAll)]
  
  res <- prcomp(msData, scale = T)
  matCP[dayInd, ] <- summary(res)$importance[3, ]
  # summary(res)
  
  pdf(paste0(PCA_folder, "/pcaVIsDrone_", the_day, ".pdf"))
  
  plot(res)
  
  op <- par(mfrow = c(2,2))
  biplot(res, choices = 1:2)
  biplot(res, choices = 3:4)
  
  factor.loadings <- cor(msData, res$x[,1:4])
  #factor.loadings
  
  theta <- 2 * pi * (0:100 / 100)
  x <- cos(theta)
  y <- sin(theta)
  plot(factor.loadings[,1:2], xlim = c(-1,1), ylim = c(-1,1), pch = " ")
  text(factor.loadings[,1:2], rownames(factor.loadings), col = "red")
  lines(x, y, col = "gray")
  abline(v = 0, h = 0)
  plot(factor.loadings[,3:4], xlim = c(-1,1), ylim = c(-1,1), pch = " ")
  text(factor.loadings[,3:4], rownames(factor.loadings), col = "red")
  lines(x, y, col = "gray")
  abline(v = 0, h = 0)
  par(op)
  
  op <- par(mfrow = c(1,2))
  my_col <- c("blue3", "darkgreen", "chocolate2", "red")
  subcondtion <- as.numeric(str_sub(dataAll$treatment, 2))
  plot(res$x[,1:2], col = my_col[subcondtion])
  legend("topright", legend = unique(dataAll$treatment), col = my_col, pch = 1)
  plot(res$x[,3:4], col = my_col[subcondtion])
  par(op)
  
  op <- par(mfrow = c(2,2))
  for (i in 1:length(unique(dataAll$treatment))) {
    # i <- 1
    each_df <- dataAll[dataAll$treatment == unique(dataAll$treatment)[i], ]
    # each_value <- each_df[, 1:6]
    dryWeightEach0 <- 2 * scale(each_df["dryWeight"], 
                                center = min(each_df["dryWeight"]), 
                                scale = max(each_df["dryWeight"]) - min(each_df["dryWeight"]))
    dryWeightEach = round(dryWeightEach0, 2)
    plot_size <- dryWeightEach + 0.3
    
    res_condition <- res$x[dataAll$treatment == unique(dataAll$treatment)[i], ]
    plot(res_condition[,1:2], col = my_col[i], cex = plot_size, pch = 19, 
         main = unique(dataAll$treatment)[i])
    #plot(res_condition[,3:4], col = my_col[i], cex = plot_size, pch = 19)
  }
  par(op)
  dev.off()
}

matCP <- round(matCP, 2)
write.csv(x = matCP, file = paste0(PCA_folder, "/cumulativeProportion.csv"))

arrayCor

###### calvulate genomic heritability #######  
correlation_folder <- paste0(PCA_folder, "/df_dry_weight_correlation")
heritability_folder <- paste0(PCA_folder, "/df_heritability")

if (!file.exists(correlation_folder)) {
  dir.create(correlation_folder)
}
if (!file.exists(heritability_folder)) {
  dir.create(heritability_folder)
}

# read the amat
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
# amat <- as.matrix(amat0)
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"

# cl <- makeCluster(15)
# registerDoParallel(cl)
foreach(dayInd = 1:length(dataWeek)) %do% {
  # dayInd <- 6
  the_day <- dataWeek[dayInd]
  dataAll0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"), 
                       header = T, row.names = 1)
  dataAll0 <- dataAll0[dataAll0$variety %in% colnames(amat0), ]
  
  # remove "plot", "flower", and each spectral
  dataAll <- dataAll0[, -c(3, 4, 7:11)]
  dataAll <- na.omit(dataAll)
  
  phenoData <- dataAll[, 3:ncol(dataAll)]
  msData <- dataAll[, 5:ncol(dataAll)]
  
  # make the mat to input the genomic heritability
  condition <- c("W1", "W2", "W3", "W4")
  h2Mat <- matrix(NA, nrow = length(condition), ncol = ncol(phenoData))
  rownames(h2Mat) <- condition
  colnames(h2Mat) <- colnames(phenoData)
  
  
  for (conditionEachInd in 1:length(condition)) {
    # conditionEachInd <- 3
    conditionEach <- condition[conditionEachInd]
    dataEachCondition <- dataAll[dataAll$treatment == conditionEach, ]
    phenoDataEachCondition <- dataEachCondition[, 3:ncol(dataEachCondition)]
    phenoDataEachConditionScaled <- scale(phenoDataEachCondition)
    rownames(phenoDataEachConditionScaled) <- dataEachCondition$variety
    amatEach <- amat0[dataEachCondition$variety, dataEachCondition$variety]
    
    # make the correlation between dryweight
    write.csv(cor(phenoDataEachConditionScaled), 
              paste0(correlation_folder, "/", the_day, "_dryWeightCor", conditionEach, ".csv"))
    
    png(paste0(correlation_folder, "/", the_day, "_Cor", conditionEach, ".png"), 
        height = 1440, width = 1440, res = 216)
    psych::pairs.panels(phenoDataEachConditionScaled)
    dev.off()
    
    phenoDataEachConditionLog <- log(phenoDataEachCondition[, c("dryWeight", "plantArea")])
    phenoDataEachConditionLogScaled <- scale(cbind(phenoDataEachConditionLog, phenoDataEachCondition[, -c(1, 2)]))
    rownames(phenoDataEachConditionLogScaled) <- dataEachCondition$variety
    write.csv(cor(phenoDataEachConditionLogScaled), 
              paste0(correlation_folder, "/", the_day, "_dryWeightCorLog", conditionEach, ".csv"))
    
    png(paste0(correlation_folder, "/", the_day, "_CorLog", conditionEach, ".png"), 
        height = 1440, width = 1440, res = 216)
    psych::pairs.panels(phenoDataEachConditionLogScaled)
    dev.off()
    
    # make amat list for each trait data
    amatZ <- design.Z(pheno.labels = rownames(phenoDataEachConditionLogScaled),
                      geno.names = rownames(amatEach))
    amatList <- list(Z = amatZ,
                     K = amatEach)
    # set each ZETA 
    ZETA <- list(amatList = amatList)
    
    for (targetTraitEachInd in 1:ncol(phenoDataEachConditionScaled)) {
      # targetTraitEachInd <- 1
      targetTraitEach <- phenoDataEachConditionLogScaled[, targetTraitEachInd]
      names(targetTraitEach) <- rownames(phenoDataEachConditionLogScaled)
      # mmfit <- mixed.solve(targetTraitEach, K = amatEach)
      mmfit <- EMM.cpp(y = targetTraitEach, ZETA = ZETA)
      h2Mat[conditionEachInd, targetTraitEachInd] <- mmfit$Vu / (mmfit$Vu + mmfit$Ve)
    }
  }
  write.csv(h2Mat, paste0(heritability_folder, "/", the_day, "_gHeritabilityEachCondition.csv"))
  
  # all condition pooled
  # calculate the narrow sense of heritability
  amatDiag <- diag(ncol(amat0))
  rownames(amatDiag) <- colnames(amatDiag) <- rownames(amat0)
  amat <- amatDiag[dataAll$variety, dataAll$variety]
  
  conditionEach <- dataAll$treatment
  treatmentMat <- model.matrix( ~ dataAll$treatment - 1)
  rownames(treatmentMat) <- dataAll$variety
  
  phenoDataScaled <- scale(phenoData)
  rownames(phenoDataScaled) <- dataAll$variety
  
  # make the correlation between dryweight
  phenoCor <- cor(phenoDataScaled)
  write.csv(phenoCor, 
            paste0(correlation_folder, "/", the_day, "_dryWeightCorAll.csv"))
  
  png(paste0(correlation_folder, "/", the_day, "_CorAllPooled.png"), 
      height = 1440, width = 1440, res = 216)
  psych::pairs.panels(phenoDataScaled)
  dev.off()
  
  phenoDataLog <- log(phenoData[, c("dryWeight", "plantArea")])
  phenoDataLogScaled <- scale(cbind(phenoDataLog, phenoData[, -c(1, 2)]))
  write.csv(cor(phenoDataLogScaled), 
            paste0(correlation_folder, "/", the_day, "_dryWeightCorAllLog.csv"))
  png(paste0(correlation_folder, "/", the_day, "_CorAllPooledLog.png"), 
      height = 1440, width = 1440, res = 216)
  psych::pairs.panels(phenoDataLogScaled)
  dev.off()
  
  # make the mat to input the narrow sense of heritability
  h2Pool <- matrix(NA, nrow = 1, ncol = ncol(phenoDataScaled))
  rownames(h2Pool) <- "All"
  colnames(h2Pool) <- colnames(phenoDataScaled)
  
  # make amat list for each trait data
  amatZ <- design.Z(pheno.labels = rownames(phenoDataScaled),
                    geno.names = rownames(amat))
  amatList <- list(Z = amatZ,
                   K = amat)
  # set each ZETA 
  ZETA <- list(amatList = amatList)
  
  for (targetTraitInd in 1:ncol(phenoDataScaled)) {
    # targetTraitInd <- 1
    X <- treatmentMat
    targetTrait <- phenoDataScaled[, targetTraitInd]
    
    # mmfit <- mixed.solve(targetTrait, K = amat, X = X)
    mmfit <- EMM.cpp(targetTrait, X = X, ZETA = ZETA, n.core = 1)
    h2Pool[, targetTraitInd] <- mmfit$Vu / (mmfit$Vu + mmfit$Ve)
  }
  write.csv(h2Pool, paste0(heritability_folder, "/", the_day, "_h2NarrowAllPooled.csv"))
}
# stopCluster(cl)

###### visualize the genomic heritability and dryWeight Correlation#######  
condition <- c("W1", "W2", "W3", "W4")

png(paste0(correlation_folder, "/dryWeightCorLog.png"), 
    height = 1440, width = 1440, res = 216)
opar <- par(mfrow = c(2, 2))
dryCorDfList <- foreach(conditionEachInd = 1:length(condition)) %do% {
  # conditionEachInd <- 1
  conditionEach <- condition[conditionEachInd]
  dryCorFileList0 <- list.files(correlation_folder, 
                                pattern = paste0("dryWeightCorLog", conditionEach), 
                                full.names = T)
  eachDayDryCorList <- lapply(dryCorFileList0, function(eachDayCorFile) {
    # eachDayCorFile <- dryCorFileList0[[1]]
    eachDayCorDf <- read.csv(eachDayCorFile, header = T, row.names = 1)
    eachDayDryCor <- eachDayCorDf["dryWeight", ]
    return(eachDayDryCor)
  })
  dryCorDf <- do.call(what = rbind, eachDayDryCorList)[, -c(1, 2)]
  rownames(dryCorDf) <- dataWeek
  corrplot(as.matrix(dryCorDf), method = "shade", tl.srt = 45, tl.cex = .8, 
           addCoef.col = "black", tl.col = "black", number.cex = .7)
  mtext(paste0("dryWeight correlation ", conditionEach)
        , at = 3, line = 3, cex = .8)
  dryCorDfEach <- cbind(conditionEach, dataWeek, dryCorDf)
  return(dryCorDfEach)
}
dev.off()
par(opar)
dryCorAllDf <- do.call(what = rbind, dryCorDfList)


png(paste0(correlation_folder, "/dryWeightCorAllLog.png"), 
    height = 1440, width = 1440, res = 216)
dryCorFileList0 <- list.files(correlation_folder, 
                              pattern = paste0("dryWeightCorAllLog.csv"), 
                              full.names = T)
eachDayDryCorList <- lapply(dryCorFileList0, function(eachDayCorFile) {
  # eachDayCorFile <- dryCorFileList0[[1]]
  eachDayCorDf <- read.csv(eachDayCorFile, header = T, row.names = 1)
  eachDayDryCor <- eachDayCorDf["dryWeight", ]
  return(eachDayDryCor)
})
dryCorDf <- do.call(what = rbind, eachDayDryCorList)[, -c(1, 2)]
rownames(dryCorDf) <- dataWeek
corrplot(as.matrix(dryCorDf), method = "shade", tl.srt = 45, tl.cex = .8, 
         addCoef.col = "black", tl.col = "black", number.cex = .7)
mtext(paste0("dryWeight correlation all log pooled")
      , at = 3, line = 3, cex = .8)
dev.off()
dryCorDfPooled <- cbind("All", dataWeek, dryCorDf)
colnames(dryCorDfPooled)[1] <-  colnames(dryCorAllDf)[1]
dryCorLogAllDf <- rbind(dryCorAllDf, dryCorDfPooled)
rownames(dryCorLogAllDf) <- NULL
# dryCorLogAllDf <- dryCorLogAllDf[, -3]
dryCorLogAllDf[, 3:ncol(dryCorLogAllDf)] <- round(dryCorLogAllDf[, 3:ncol(dryCorLogAllDf)], 2)
write.csv(dryCorLogAllDf, paste0(correlation_folder, "/dryWeightCorLogAllCondition.csv"))



# visualize the genomic heritability
gHerFileList <- list.files(heritability_folder, 
                           pattern = "gHeritability", 
                           full.names = T)
png(paste0(heritability_folder, "/genomicHeritability.png"), 
    height = 1440, width = 1440, res = 216)
opar <- par(mfrow = c(2, 2))
herDfList <- foreach(conditionEachInd = 1:length(condition)) %do% {
  # conditionEachInd <- 1
  conditionEach <- condition[conditionEachInd]
  eachDayHerList <- lapply(gHerFileList, function(eachDayHerFile) {
    # eachDayHerFile <- gHerFileList[[1]]
    eachDayHerDf <- read.csv(eachDayHerFile, header = T, row.names = 1)
    eachDayHer <- eachDayHerDf[conditionEach, ]
    return(eachDayHer)
  })
  herDf <- do.call(what = rbind, eachDayHerList)[, -c(1, 2)]
  rownames(herDf) <- dataWeek
  corrplot(as.matrix(herDf), method = "shade", tl.srt = 45, tl.cex = .8, 
           addCoef.col = "black", tl.col = "black", number.cex = .7)
  mtext(paste0("genomic heritability ", conditionEach)
        , at = 3, line = 3, cex = .8)
  herDfEach <- cbind(conditionEach, dataWeek, herDf)
  return(herDfEach)
}
dev.off()
par(opar)
herAllDf <- do.call(what = rbind, herDfList)
rownames(herAllDf) <- NULL
# herAllDf <- herAllDf[, -c(3, 4)]
herAllDf[, 3:ncol(herAllDf)] <- round(herAllDf[, 3:ncol(herAllDf)], 2)
write.csv(herAllDf, paste0(heritability_folder, "/genomicHeritabilityLogAllCondition.csv"))


# visualize the narrow sense of heritability
h2NarrowFileList <- list.files(heritability_folder, 
                               pattern = "h2Narrow", 
                               full.names = T)
png(paste0(heritability_folder, "/narrowHeritability.png"), 
    height = 1440, width = 1440, res = 216)
eachDayHerList <- lapply(h2NarrowFileList, function(eachDayHerFile) {
  # eachDayHerFile <- h2NarrowFileList[[1]]
  eachDayHerDf <- read.csv(eachDayHerFile, header = T, row.names = 1)
  return(eachDayHerDf)
})
herDf <- do.call(what = rbind, eachDayHerList)
rownames(herDf) <- dataWeek
corrplot(as.matrix(herDf), method = "shade", tl.srt = 45, tl.cex = .8, 
         addCoef.col = "black", tl.col = "black", number.cex = .7)
mtext(paste0("narrow heritability"))
dev.off()
par(opar)

####### make box plot of dryWeight in each condition #########
phenoAll0 <- read.csv(paste0(phenoFolder, "/week1_medianPheno.csv"), 
                      header = T, row.names = 1)
phenoAll <- na.omit(phenoAll0)
phenoAll <- phenoAll[phenoAll$variety %in% colnames(amat0), ]

png(paste0(PCA_folder, "/boxplotBiomass.png"), 
    height = 1440, width = 1440, res = 216)
# boxplot(dryWeight ~ treatment, data = phenoAll, 
#         ylab = "dry weight (g)", main = "boxplot of dry weight", 
#         col = c("blue3", "lightskyblue1", "darkgoldenrod1", "chocolate2"))

g <- ggplot(phenoAll, aes(x = treatment, y = dryWeight, fill = treatment)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("blue3", "lightskyblue1", "darkgoldenrod1", "chocolate2")) + 
  labs(title = "boxplot of dry weight", y = "dry weight (g)")
print(g)
dev.off()

# t.test
phenoAll$treatment <- as.factor(phenoAll$treatment)
anova(lm(phenoAll$dryWeight~phenoAll$treatment))

t.test(phenoAll[phenoAll$treatment == "W1", "dryWeight"], 
       phenoAll[phenoAll$treatment == "W2", "dryWeight"])$p.value * 6
t.test(phenoAll[phenoAll$treatment == "W1", "dryWeight"], 
       phenoAll[phenoAll$treatment == "W3", "dryWeight"])$p.value * 6
t.test(phenoAll[phenoAll$treatment == "W1", "dryWeight"], 
       phenoAll[phenoAll$treatment == "W4", "dryWeight"])$p.value * 6
t.test(phenoAll[phenoAll$treatment == "W2", "dryWeight"], 
       phenoAll[phenoAll$treatment == "W3", "dryWeight"])$p.value * 6
t.test(phenoAll[phenoAll$treatment == "W2", "dryWeight"], 
       phenoAll[phenoAll$treatment == "W4", "dryWeight"])$p.value * 6
t.test(phenoAll[phenoAll$treatment == "W3", "dryWeight"], 
       phenoAll[phenoAll$treatment == "W4", "dryWeight"])$p.value * 6


# calculate the MEAN and SE
dryWeightMean <- tapply(X = phenoAll$dryWeight, INDEX = phenoAll$treatment, mean)
dryWeightSd <- tapply(X = phenoAll$dryWeight, INDEX = phenoAll$treatment, sd)
dryWeightLength <- tapply(X = phenoAll$dryWeight, INDEX = phenoAll$treatment, length)
dryWeightSe <- dryWeightSd / sqrt(dryWeightLength)
dryWeightData <- rbind(dryWeightMean, dryWeightSd, dryWeightSe)
rownames(dryWeightData) <- c("mean", "SD", "SE")
dryWeightData <- round(dryWeightData, 2)
write.csv(dryWeightData, paste0(correlation_folder, "/dryWeightSummary.csv"))
