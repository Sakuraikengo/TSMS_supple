# machine <- "drone"
the_year <- "2019"
data_day <- c("0802", "0810", "0817", "0824", "0831", "0903")
source("scripts/1.functionCode.R")

options(stringsAsFactors = FALSE)

library(readr)
library(stringr)
library(psych)
library(ggplot2)
# library(gaston)

# make the folder
resultFolder <- "2019_M100_Xacti4eyeCamera_Images/result/9.0.0.ANOVAdryweight"
if (!file.exists(resultFolder)) {
  dir.create(resultFolder)
}

# fill NA irregular plot
irr_plot <- c("W1-036", "W1-073", "W1-085", "W1-115", "W1-129",
              "W2-005", "W2-038", "W2-099", "W2-149", "W2-051",
              "W3-014", "W3-036", "W3-080", "W3-093", "W3-114", 
              "W4-004", "W4-071", "W4-074", "W4-075", "W4-150", "W4-182", "W4-184", 
              "W4-009", "W4-143", "W4-034", "W4-042", "W4-109")


# read the genome data
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"

# read the dry weight
plantPhenotype0 <- read.csv("2019_M100_Xacti4eyeCamera_Images/field_data/2019_Tottori_Jul_ShootPhenotype.csv", 
                            header = TRUE, row.names = 1, stringsAsFactors = T)
start_W <- match("W1", plantPhenotype0$block)
plantPhenotype <- plantPhenotype0[start_W:dim(plantPhenotype0)[1], ]
plantPhenotype <- plantPhenotype[plantPhenotype$plot <= 200, ]
plantPhenotype$year <- as.factor(plantPhenotype$year)
plantPhenotype$plot <- as.factor(plantPhenotype$plot)
plantPhenotype$varietyID <- as.factor(plantPhenotype$varietyID)
plantPhenotype$ind <- as.factor(plantPhenotype$ind)
plotNumber <- paste0(plantPhenotype$block, "-", formatC(plantPhenotype$plot, width = 3, flag = "0"))
plantPhenotype <- cbind(plantPhenotype, plotNumber)

plantPhenotype[plantPhenotype$plotNumber %in% irr_plot, 9:(ncol(plantPhenotype) - 1)] <- NA
plantPhenotype[!(plantPhenotype$variety %in% colnames(amat0)), 9:(ncol(plantPhenotype) - 1)] <- NA

plantWeight <- na.omit(plantPhenotype[, c(1:8, 10)])
plantWeight$variety <- as.character(plantWeight$variety)
plantWeight$variety <- as.factor(plantWeight$variety)
# summary(plantPhenotype)

# calculate the anova in each treatment (the effect of accessions)
anova(lm(DryWeight_Shoot_g ~ block + variety, 
         data = plantWeight))
# anova(lm(DryWeight_Shoot_g ~ plot,
#          data = plantPhenotype[plantPhenotype$block == "W1", ]))
# anova(lm(DryWeight_Shoot_g ~ plot, 
#          data = plantPhenotype[plantPhenotype$block == "W2", ]))
# anova(lm(DryWeight_Shoot_g ~ plot, 
#          data = plantPhenotype[plantPhenotype$block == "W3", ]))
# anova(lm(DryWeight_Shoot_g ~ plot, 
#          data = plantPhenotype[plantPhenotype$block == "W4", ]))

# dryLeaves <- aggregate(DryWeight_Leaves_g~plotNumber, data = plantPhenotype, FUN = mean)
dryShoot <- aggregate(DryWeight_Shoot_g~plotNumber, data = plantPhenotype, FUN = mean)

# dryLeavesInd <- match(x = unique(plotNumber), table = dryLeaves$plotNumber)
dryShootInd <- match(x = unique(plotNumber), table = dryShoot$plotNumber)
# dryWeight0 <- dryLeaves[dryLeavesInd, "DryWeight_Leaves_g"] + dryShoot[dryShootInd, "DryWeight_Shoot_g"]
dryWeight0 <- dryShoot[dryShootInd, "DryWeight_Shoot_g"]

variety <- plantPhenotype$variety[seq(1,nrow(plantPhenotype), 4)]
dryWeight <- data.frame(variety = variety, 
                        block = str_sub(unique(plotNumber), 1, 2), 
                        dryWeight = dryWeight0)
rownames(dryWeight) <- unique(plotNumber)

# calculate the mean and standard deviation
se <- function(x) sd(x) / sqrt(length(x))
test <- na.omit(dryWeight)
print(tapply(test$dryWeight, test$block, se))
print(tapply(test$dryWeight, test$block, mean))

# read the flowering date
flowerDay0 <- read.csv("2019_M100_Xacti4eyeCamera_Images/field_data/2019_Tottori_Jul_FloweringDate.csv", 
                       header = T)
flowerDay <- flowerDay0[, c("variety", "block", "FloweringDate")]
flowerDate <- FlowerDataAsDate(flowerDate = flowerDay)
rownames(flowerDate) <- rownames(dryWeight)
dayInd <- as.Date(c("2019-08-02", "2019-08-10", "2019-08-17", 
                    "2019-08-24", "2019-08-31", "2019-09-03"))

# calculate the correlation between DTF and dry weight
dfDD <- data.frame(dryWeight, 
                   DTF = as.numeric(flowerDate$FloweringDate))
cor(dfDD[1:200, 3:4], use = "complete.obs")[1,2]
cor(dfDD[201:400, 3:4], use = "complete.obs")[1,2]
cor(dfDD[401:600, 3:4], use = "complete.obs")[1,2]
cor(dfDD[601:800, 3:4], use = "complete.obs")[1,2]

# save the DTF
write.csv(dfDD, paste0(resultFolder, "/DTF.csv"))

flowerDateList <- lapply(dayInd, function(eachDayInd) {
  # eachDayInd <- dayInd[1]
  flowerOrNot <- flowerDate[, 3] - eachDayInd
  # 0 means flowering already, 1 means not flowering, NA is also not flowering
  flowerOrNot[flowerOrNot > 0] <- 1
  flowerOrNot[flowerOrNot <= 0] <- 0
  flowerOrNot[is.na(flowerOrNot)] <- 1
  return(flowerOrNot)
})

flowerMat <- do.call(what = cbind, args = flowerDateList)
pheno <- data.frame(dryWeight, flowerMat)
colnames(pheno) <- c("variety", "block", "dryWeight", 
                     "Week1", "Week2", "Week3", "Week4", "Week5", "Week6")
pheno$block <- as.factor(pheno$block)
pheno$Week1 <- as.factor(pheno$Week1)
pheno$Week2 <- as.factor(pheno$Week2)
pheno$Week3 <- as.factor(pheno$Week3)
pheno$Week4 <- as.factor(pheno$Week4)
pheno$Week5 <- as.factor(pheno$Week5)
pheno$Week6 <- as.factor(pheno$Week6)

summary(pheno)
week <- paste0("Week", 2:6)
# calculate the anova (the effect of flowering date)
# anova(lm(dryWeight ~ Week1, data = pheno))
anovaList <- lapply(week, function(eachWeek) {
  phenoEach <- na.omit(pheno)
  colnames(phenoEach)[colnames(phenoEach) == eachWeek] <- "flower"
  # f0 <- as.matrix(anova(lm(dryWeight ~ block + flower, data = phenoEach)))
  f1 <- as.matrix(anova(lm(dryWeight ~ flower,
                           data = phenoEach[phenoEach$block == "W1", ])))
  f2 <- as.matrix(anova(lm(dryWeight ~ flower,
                           data = phenoEach[phenoEach$block == "W2", ])))
  f3 <- as.matrix(anova(lm(dryWeight ~ flower,
                           data = phenoEach[phenoEach$block == "W3", ])))
  f4 <- as.matrix(anova(lm(dryWeight ~ flower,
                           data = phenoEach[phenoEach$block == "W4", ])))
  f <- rbind(f1, f2, f3, f4)
  rownames(f) <- rep(c(eachWeek, "Residuals"), 4)
  # rownames(f0) <- c("Treatment", eachWeek, "Residuals")
  return(f)
})
anovaMat <- do.call(rbind, anovaList)
write.csv(anovaMat, paste0(resultFolder, "/flowerANOVA.csv"))

##### calculating the effect of GxE
phenoDataAll0 <- na.omit(pheno)
phenoDataAll0 <- phenoDataAll0[phenoDataAll0$variety %in% colnames(amat0), ]

# take log
phenoDataAll0$dryWeight <- log(phenoDataAll0$dryWeight)

# meke treatment kernel
treatmentMat <- model.matrix( ~ phenoDataAll0$block - 1)
rownames(treatmentMat) <- rownames(phenoDataAll0)
kernelTre <- tcrossprod(treatmentMat)
# image(kernelTre)

# dryweight
dryWeight <- phenoDataAll0[, "dryWeight"]
names(dryWeight) <- phenoDataAll0$variety

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

# make G * treatment list
amatTreZ <- design.Z(pheno.labels = names(dryWeight),
                     geno.names = rownames(kernelAmatTre))

amatTreList <- list(Z = amatTreZ,
                    K = kernelAmatTre)

ZETA_GT <- list(amat = amatList,
                amatTre = amatTreList)
# make the fixed effect
treatmentX0 <- treatmentMat

for (i in 1:length(flowerDateList)) {
  # i <- 1
  df <- phenoDataAll0[, c(1:3, i + 3)]
  colnames(df)[ncol(df)] <- "flower"
  
  # make the fixed effect
  flowerDate <- df$flower
  
  if (any(flowerDate == 0)) {
    flowerX0 <- cbind(1, flowerDate)
    rownames(flowerX0) <- names(dryWeight)
    multiX0 <- cbind(treatmentX0, flowerDate)
  } else {
    flowerX0 <- NULL
    multiX0 <- treatmentX0
  }
  rownames(multiX0) <- names(dryWeight)
  
  ZETA <- ZETA_GT
  X0 <- multiX0
  nameInd <- "G+E+GEI"
  
  resEM3 <- EM3.cpp(y = dryWeight,
                    X0 = X0, n.core = 1,
                    ZETA = ZETA)
  str(resEM3)
}
