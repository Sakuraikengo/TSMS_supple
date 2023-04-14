# machine <- "drone"
the_year <- "2019"
dataWeek <- c("week1", "week2", "week3", "week4", "week5", "week6")

#' # 0. Read packages
library(rrBLUP)
library(MASS)
library(tidyr)
library(MTM)
library(corrplot)
library(stringr)
library(date)
library(ggplot2)
library(ggsci)
library(doParallel)
source("scripts/1.functionCode.R")
source("scripts/1.MTM_2.R")
options(stringsAsFactors = FALSE)

# make the folder
mtmFolder <- "2019_M100_Xacti4eyeCamera_Images/result/3.0.2.MTMFolderOnlyNDVI"
# mtmFolder <- "C:/Users/biometrics/Desktop/3.0.2.test"
if (!file.exists(mtmFolder)) {
  dir.create(mtmFolder)
}

# data folder
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
foreach(dayInd = 1:length(dataWeek), .packages = c("MTM", "corrplot")) %dopar% {
  # the_day <- dataWeek[[4]]
  the_day <- dataWeek[[dayInd]]
  
  dataFrame0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"),
                         header = TRUE, row.names = 1)
  df <- dataFrame0[, c("variety", "treatment", "flower", "dryWeight", "NDVI")]
  df <- na.omit(df)
  
  
  n_treatment <- length(unique(df$treatment))
  name_treatment <- c("W1", "W2", "W3", "W4")
  name_trait <- colnames(df[, c(4:ncol(df))])
  n_trait <- length(name_trait)
  
  phenodat <- df
  colnames(phenodat)[c(1, 2)] <- c("line", "treatment")
  
  
  #' # 3. Estimate genotypic values with MTM
  #' Estimate genotypic values for each condition,
  #' ohterwise the variance-covariance matrix to estimate will be too large
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
    # estimation
    MTM_output <- capture.output(
      res_MTM <- MTM_2(XF = MTM_XF, Y = MTM_Y, K = MTM_K, resCov = MTM_resCov,
                       nIter = 120000, burnIn = 20000, thin = 200, 
                       saveAt = paste0(new_folder, "/", treatment_i, '_'))
    )
    
    write.table(res_MTM$YHat, paste0(new_folder, "/", treatment_i, "_MTM_predictions.txt"), quote = F)
    write.table(res_MTM$resCov$R, paste0(new_folder, "/", treatment_i, "_MTM_resCov.txt"), quote = F)
    write.table(res_MTM$K[[1]]$G, paste0(new_folder, "/", treatment_i, "_MTM_G.txt"), quote = F)
    write.table(res_MTM$K[[1]]$U, paste0(new_folder, "/", treatment_i, "_MTM_U.txt"), quote = F)
    
    # diagnosis (log-likelihood)
    MTM_logLik <- read.table(paste0(new_folder, "/", treatment_i, '_logLik.dat'))
    png(paste0(new_folder, "/", treatment_i, "_logLik.png"))
    plot(MTM_logLik[[1]], cex = 0.5,
         main = paste('Log-likelihood, ', treatment_i, the_day),
         ylab = 'Log-likelihood')
    dev.off()
    
    # diagnosis (mu)
    MTM_mu <- readLines(paste0(new_folder, "/", treatment_i, '_mu.dat'))
    MTM_mu_1 <- sapply(MTM_mu, strsplit, split = ' ')
    names(MTM_mu_1) <- NULL
    MTM_mu_2 <- sapply(MTM_mu_1, function(x) {
      length_max <- max(sapply(MTM_mu_1, length))
      y <- rep(NA, length_max)
      y[1:length(x)] <- x
      as.numeric(y)
    })
    MTM_mu_3 <- t(as.data.frame(MTM_mu_2))
    png(paste0(new_folder, "/", treatment_i, "_mu.png"))
    matplot(MTM_mu_3, cex = 0.5, pch = 1:ncol(MTM_mu_3),
            main = paste0('mu, ', treatment_i, '_', the_year, the_day),
            ylab = 'mu')
    dev.off()
  }
}
stopCluster(cl)

#### make corplot from "genetic covariance matrix" "residual covariance matrix" ######
#### make ggplot for the change of "obs_cor", "g_cor", "res_cor" #####
band_name <- c("dryWeight", "NDVI")
condition <- c("W1", "W2", "W3", "W4")  


# make array to input the gCor between dryWeight and VIs
gCorDryWeightArray <- array(NA, dim = c(length(dataWeek), length(band_name), length(condition)))
dimnames(gCorDryWeightArray) <- list(dataWeek, band_name, condition)

# save the corplot as png
treatment <- c("C", "W5", "W10", "D")

for (the_day in dataWeek) {
  # the_day <- dataWeek[4]
  
  # set working directory
  data_folder <- paste0(mtmFolder, "/", the_year, the_day)
  
  # read the cov files
  g_cov_files <- list.files(data_folder, pattern = "MTM_G")
  res_cov_files <- list.files(data_folder, pattern = "MTM_resCov")
  
  g_cor_list <- list()
  res_cor_list <- list()
  
  png(paste0(mtmFolder, "/", the_day, "geneticCor.png"), 
      height = 1440, width = 1440, res = 216)
  opar <- par(mfrow = c(2, 2))
  
  
  for (file_index in 1:length(g_cov_files)) {
    # file_index <- 1
    treatmentEach <- treatment[file_index]
    g_cov_read <- g_cov_files[file_index]
    
    # make cor plot "g"
    g_cov_path <- paste0(data_folder, "/", g_cov_read)
    indexAll <- c("dryWeight", "NDVI")
    cor_value_g0 <- MakeCorMat(g_cov_path, indexAll)
    cor_value_g <- cor_value_g0
    
    corrplot(cor_value_g, method = "shade", tl.srt = 45, addCoef.col = "black", 
             tl.col = "black", number.cex = .7)
    mtext(paste0(str_to_title(the_day), " genetic correlation ", treatmentEach)
          , at = 3, line = 2.5, cex = 1)
    
    gCorDryWeightArray[the_day, , file_index] <- cor_value_g["dryWeight", ]
  } 
  dev.off()
  par(opar)
}

# make the day change of the genetic correlation between biomass and VIs
eachCorArray <- gCorDryWeightArray[, -1, ]
rownames(eachCorArray) <- str_replace(rownames(eachCorArray), 
                                      pattern = "w", replacement = "W")
colnames(eachCorArray) <- c("C", "W5", "W10", "D")
png(paste0(mtmFolder, "/genomeCorDryWeight.png"), 
    height = 1440, width = 1440, res = 216)
write.csv(eachCorArray, 
          paste0(mtmFolder, "/genomeCorDryWeight_NDVI.csv"))
corrplot(eachCorArray, method = "shade", tl.srt = 45, 
         addCoef.col = "black", tl.col = "black")
# mtext(paste0("genomeCorDryWeight NDVI")
#       , at = 3, line = 3, cex = .8)

dev.off()

# from week4 to week5 what % increase
mean(gCorDryWeightArray[5, 2:6, 1:3] / gCorDryWeightArray[4, 2:6, 1:3]) 
