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
options(stringsAsFactors = FALSE)

# make the folder
mtmFolder <- "2019_M100_Xacti4eyeCamera_Images/result/3.0.0.MTMFolder"
if (!file.exists(mtmFolder)) {
  dir.create(mtmFolder)
}

# data folder
phenoFolder <- "2019_M100_Xacti4eyeCamera_Images/result/2.1.1.dataFolder"

# read the amat
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"


#### calculate the MTM ######
cl <- makeCluster(15)
registerDoParallel(cl)
foreach(dayInd = 1:length(dataWeek), .packages = c("MTM", "corrplot")) %dopar% {
  # dayInd <- 1
  the_day <- dataWeek[dayInd]
  
  # read the phenotypic data
  dataFrame0 <- read.csv(paste0(phenoFolder, "/", the_day, "_medianPheno.csv"),
                         header = TRUE, row.names = 1)
  df <- dataFrame0[, c("variety", "treatment", "dryWeight", "plantArea", 
                       "GRVI",
                       "NDVI", 
                       "NDRE", 
                       "NDI", 
                       "RTVI")]
  df <- na.omit(df)
  
  # set the treatments and traits
  n_treatment <- length(unique(df$treatment))
  name_treatment <- c("W1", "W2", "W3", "W4")
  name_trait <- colnames(df[, c(3:ncol(df))])
  n_trait <- length(name_trait)
  
  phenodat <- df
  colnames(phenodat)[c(1, 2)] <- c("line", "treatment")
  
  
  # 3. Estimate genotypic values with MTM
  # Estimate genotypic values for each condition,
  new_folder <- paste0(mtmFolder, "/", the_year, the_day)
  if (!file.exists(new_folder)) {
    dir.create(new_folder)
  }
  for (treatment_i in name_treatment) {
    # treatment_i <- name_treatment[3]
    
    # extract the data of the target treatment
    phenodat_i <- phenodat[phenodat$treatment == treatment_i, ]
    phenodat_i <- na.omit(phenodat_i)
    
    
    phenodat_i <- phenodat_i[phenodat_i$line %in% colnames(amat0), ]
    amat <- amat0[as.vector(phenodat_i$line), as.vector(phenodat_i$line)]
    name_line <- rownames(amat)
    
    # Y
    MTM_Y <- as.matrix(phenodat_i[, -c(1:2)])
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
    
    # estimation
    MTM_output <- capture.output(
      res_MTM <- MTM(Y = MTM_Y, K = MTM_K, resCov = MTM_resCov,
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
band_name <- c("dryWeight", 
               "GRVI",
               "NDVI", 
               "NDRE", 
               "NDI", 
               "RTVI")
condition <- c("W1", "W2", "W3", "W4")  


# make array to input the gCor between dryWeight and VIs
gCorDryWeightArray <- array(NA, dim = c(length(dataWeek), length(band_name), length(condition)))
dimnames(gCorDryWeightArray) <- list(dataWeek, band_name, condition)

# save the corplot as png
treatment <- c("C", "W5", "W10", "D")

for (the_day in dataWeek) {
  # the_day <- dataWeek[1]
  
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
    indexAll <- c("dryWeight", "plantArea", "GRVI", "NDVI", "NDRE", "NDI", "RTVI")
    cor_value_g0 <- MakeCorMat(g_cov_path, indexAll)
    cor_value_g <- cor_value_g0[-2, -2]
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
png(paste0(mtmFolder, "/genomeCorDryWeight.png"), 
    height = 1440, width = 1440, res = 216)
opar <- par(mfrow = c(2, 2))

for (conditionEach in condition) {
  # conditionEach <- condition[1]
  write.csv(eachCorArray[, , conditionEach], 
            paste0(mtmFolder, "/genomeCorDryWeight_", conditionEach, ".csv"))
  corrplot(eachCorArray[, , conditionEach], method = "shade", tl.srt = 45, tl.cex = .8, 
           addCoef.col = "black", tl.col = "black", number.cex = .7)
  mtext(paste0("genomeCorDryWeight ", conditionEach)
        , at = 3, line = 3, cex = .8)
}
dev.off()
par(opar)
