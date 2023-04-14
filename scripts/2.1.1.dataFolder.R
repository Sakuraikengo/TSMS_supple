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
resultFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"
if (!file.exists(resultFolder)) {
  dir.create(resultFolder)
}

# fill NA irregular plot
irr_plot <- c("W1-036", "W1-073", "W1-085", "W1-115", "W1-129",
              "W2-005", "W2-038", "W2-099", "W2-149", "W2-051",
              "W3-014", "W3-036", "W3-080", "W3-093", "W3-114", 
              "W4-004", "W4-071", "W4-074", "W4-075", "W4-150", "W4-182", "W4-184", 
              "W4-009", "W4-143", "W4-034", "W4-042", "W4-109")

# read the plant area
plantAreaAll0 <- read.csv("2019_M100_Xacti4eyeCamera_Images/plantArea/allPlantArea.csv", 
                          header = T, row.names = 1)
plantAreaAll <- plantAreaAll0[order(rownames(plantAreaAll0)), ]


# read the dry weight
plantPhenotype <- read.csv("2019_M100_Xacti4eyeCamera_Images/field_data/2019_Tottori_Jul_ShootPhenotype.csv", 
                           header = TRUE, row.names = 1)
start_W <- match("W1", plantPhenotype$block)
plantPhenotype <- plantPhenotype[start_W:dim(plantPhenotype)[1], ]

plotNumber <- paste0(plantPhenotype$block, "-", formatC(plantPhenotype$plot, width = 3, flag = "0"))
plantPhenotype <- cbind(plantPhenotype, plotNumber)

# dryLeaves <- aggregate(DryWeight_Leaves_g~plotNumber, data = plantPhenotype, FUN = mean)
dryShoot <- aggregate(DryWeight_Shoot_g~plotNumber, data = plantPhenotype, FUN = mean)
# dryLeavesInd <- match(x = rownames(plantAreaAll), table = dryLeaves$plotNumber)
dryShootInd <- match(x = rownames(plantAreaAll), table = dryShoot$plotNumber)
dryWeight0 <- dryShoot[dryShootInd, "DryWeight_Shoot_g"]

dryWeight <- matrix(NA, nrow = nrow(plantAreaAll), ncol = 1)
rownames(dryWeight) <- rownames(plantAreaAll)
dryWeight[, 1] <- dryWeight0


medianCsvList <- list.files("2019_M100_Xacti4eyeCamera_Images/MSdata",
                            pattern = "median", full.names = T)
eachDayDfList <- lapply(medianCsvList, function(eachdayCsv) {
  # eachdayCsv <- medianCsvList[[1]]
  
  # read the plant area each day
  theDay0 <- str_split(eachdayCsv, pattern = "/")[[1]][3]
  theDay <- str_sub(theDay0, 1, 4)
  plantAreaCol <- paste0("plantArea_", theDay)
  plantAreaEachDay <- plantAreaAll[, plantAreaCol]
  
  # read the each day MS data
  eachDayDf0 <- read.csv(eachdayCsv, header = T, row.names = 1)
  
  # remove the 300s plot
  plotInd <- as.numeric(str_sub(rownames(eachDayDf0), start = 4, end = 6))
  eachDayDf <- eachDayDf0[plotInd <= 200, ]
  spectral475 <- eachDayDf$X475550850_475
  spectral550 <- (eachDayDf$X475550850_550 + eachDayDf$X550660850_550) / 2
  spectral660 <- eachDayDf$X550660850_660
  spectral725 <- eachDayDf$X725
  spectral850 <- (eachDayDf$X475550850_850 + eachDayDf$X550660850_850) / 2
  eachDayMs <- cbind(spectral475, spectral550, spectral660, spectral725, spectral850,
                     GRVI = eachDayDf$GRVI,
                     NDVI = eachDayDf$NDVI, 
                     NDRE = eachDayDf$NDRE, 
                     NDI = eachDayDf$NDI, 
                     RTVI = eachDayDf$RTVI)
  rownames(eachDayMs) <- rownames(eachDayDf)
  
  eachDayPhenoDf <- cbind(dryWeight = dryWeight[, 1], plantArea = plantAreaEachDay, eachDayMs)
  eachDayPhenoDf[rownames(eachDayPhenoDf) %in% irr_plot, ] <- NA
  return(eachDayPhenoDf)
})

# read the flowering date
flowerDay0 <- read.csv("2019_M100_Xacti4eyeCamera_Images/field_data/2019_Tottori_Jul_FloweringDate.csv", 
                       header = T)
flowerDay <- flowerDay0[, c("variety", "block", "FloweringDate")]
flowerDate <- FlowerDataAsDate(flowerDate = flowerDay)
rownames(flowerDate) <- rownames(plantAreaAll)
dayInd <- as.Date(c("2019-08-02", "2019-08-10", "2019-08-17", 
                    "2019-08-24", "2019-08-31", "2019-09-03", 
                    "2019-09-11"))

flowerDateList <- lapply(dayInd, function(eachDayInd) {
  # eachDayInd <- dayInd[1]
  flowerOrNot <- flowerDate[, 3] - eachDayInd
  # 0 means flowering already, 1 means not flowering, NA is also not flowering
  flowerOrNot[flowerOrNot > 0] <- 1
  flowerOrNot[flowerOrNot <= 0] <- 0
  flowerOrNot[is.na(flowerOrNot)] <- 1
  return(flowerOrNot)
})

data_day <- c("0802", "0810", "0817", "0824", "0831", "0903")
day1Df <- cbind(flower = flowerDateList[[1]], 
                eachDayDfList[[1]])
day2Df <- cbind(flower = flowerDateList[[2]], 
                eachDayDfList[[2]])
day3Df <- cbind(flower = flowerDateList[[3]], 
                eachDayDfList[[3]])
day4Df <- cbind(flower = flowerDateList[[4]], 
                eachDayDfList[[4]])
day5Df <- cbind(flower = flowerDateList[[5]], 
                eachDayDfList[[5]])
day6Df <- cbind(flower = flowerDateList[[6]], 
                eachDayDfList[[6]])


dataWeek <- c("week1", "week2", "week3", "week4", "week5", "week6")
allDayList <- list(list(df = day1Df, name = "week1"), 
                   list(df = day2Df, name = "week2"), 
                   list(df = day3Df, name = "week3"), 
                   list(df = day4Df, name = "week4"), 
                   list(df = day5Df, name = "week5"), 
                   list(df = day6Df, name = "week6"))

# save the phenotype
lapply(allDayList, function(eachDayList) {
  # eachDayList <- allDayList[[1]]
  eachDayDfAll <- cbind(plantAreaAll[, 1:3], eachDayList$df)
  write.csv(eachDayDfAll, 
            file = paste0(resultFolder, "/", eachDayList$name, "_medianPheno.csv"))
  
})

# percentage of flowering or not in each measurement day and each treatment
flowerRateList <- lapply(allDayList, function(eachDayList) {
  # eachDayList <- allDayList[[1]]
  eachDayDfAll <- cbind(plantAreaAll[, 1:3], eachDayList$df)
  eachDayDfAll <- na.omit(eachDayDfAll)
  flowerRateEach <- tapply(eachDayDfAll$flower, eachDayDfAll$treatment, function(eachTre) {
    (length(eachTre) - sum(eachTre)) / length(eachTre)
  })
  flowerRateAll <- (length(eachDayDfAll$flower) - sum(eachDayDfAll$flower)) / length(eachDayDfAll$flower)
  flowerRate <- c(flowerRateAll, flowerRateEach)
  return(flowerRate)
})
flowerRateMat0 <- do.call(what = rbind, flowerRateList)

# add the destructive measurement
lastFlower <- flowerDateList[[length(flowerDateList)]]
ind <- rep(c("C", "W5", "W10", "D"), each = 200)
flowerRateEach <- tapply(lastFlower, ind, function(eachTre) {
  (length(eachTre) - sum(eachTre)) / length(eachTre)
})
flowerRateAll <- (length(lastFlower) - sum(lastFlower)) / length(lastFlower)
flowerRate <- c(flowerRateAll, flowerRateEach)

flowerRateMat <- rbind(flowerRateMat0, flowerRate)
flowerRateMat <- round(flowerRateMat, 2)
rownames(flowerRateMat) <- c(dataWeek, "Destructive")
colnames(flowerRateMat) <- c("All", "C", "W5", "W10", "D")
write.csv(flowerRateMat, paste0(resultFolder, "/floweringRate.csv"))

########### visualize flowering time###########
sowingDay <- as.Date("2019-07-10")
measuringDay <- as.Date(c("2019-08-02", "2019-08-10", "2019-08-17", 
                          "2019-08-24", "2019-08-31", "2019-09-03"))

DTF <- flowerDate$FloweringDate - sowingDay
flowerDf <- cbind(flowerDate, DTF)
flowerDf <- na.omit(flowerDf)
g <- ggplot(flowerDf, aes(x = FloweringDate)) + 
  facet_wrap(~block) + 
  geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = measuringDay, linetype = "dashed") + 
  scale_x_date(date_breaks = "1 months") + 
  xlab("Flowering Day") + 
  ylab("Count") + 
  theme(plot.title = element_text(size = 16, face = "bold"), 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12,face = "bold"), 
        legend.text = element_text(size = 16))

png(paste0(resultFolder, "/floweringDay.png"), 
    height = 1440, width = 1440, res = 216)
print(g)
dev.off()

