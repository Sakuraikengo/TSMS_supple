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
mtmFolder <- "2019_M100_Xacti4eyeCamera_Images/result/9.0.1.MTMFolderOnlyDTF/"
# mtmFolder <- "C:/Users/biometrics/Desktop/9.0.1.test"
if (!dir.exists(mtmFolder)) {
  dir.create(mtmFolder)
}

# fill NA irregular plot
irr_plot <- c("W1-036", "W1-073", "W1-085", "W1-115", "W1-129",
              "W2-005", "W2-038", "W2-099", "W2-149", "W2-051",
              "W3-014", "W3-036", "W3-080", "W3-093", "W3-114", 
              "W4-004", "W4-071", "W4-074", "W4-075", "W4-150", "W4-182", "W4-184", 
              "W4-009", "W4-143", "W4-034", "W4-042", "W4-109")

# read the DTF file
phenoFile <- "2019_M100_Xacti4eyeCamera_Images/result/9.0.0.ANOVAdryweight/DTF.csv"
phenoData0 <- read.csv(phenoFile, header = T, row.names = 1)

# DTF == NA means that this plot did not flower, so set the NA as final measurement day
phenoData0[is.na(phenoData0$DTF), "DTF"] <- max(phenoData0$DTF, na.rm = T)

# remove the irr plot
df <- phenoData0[!(rownames(phenoData0) %in% irr_plot), ]
df <- na.omit(df)

# read the amat
amat0 <- as.matrix(read.csv('genome/amat173583SNP.csv', 
                            row.names = 1, header = T))
# amat <- as.matrix(amat0)
rownames(amat0)[rownames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "HOUJAKU_KUWAZU"] <- "Houjaku Kuwazu"
colnames(amat0)[colnames(amat0) == "X5002T"] <- "5002T"

# set the data information
name_treatment <- c("W1", "W2", "W3", "W4")
name_trait <- colnames(df[, c(3:ncol(df))])
n_trait <- length(name_trait)

phenodat <- df
colnames(phenodat)[c(1, 2)] <- c("line", "treatment")


#' # 3. Estimate genotypic values with MTM
#' Estimate genotypic values for each condition,
#' ohterwise the variance-covariance matrix to estimate will be too large
cl <- makeCluster(4)
registerDoParallel(cl)

foreach(treatment_i = name_treatment, .packages = c("MTM", "corrplot")) %dopar% {
  # treatment_i <- name_treatment[1]
  
  # choose data of the target treatment
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
    res_MTM <- MTM_2(XF = NULL, Y = MTM_Y, K = MTM_K, resCov = MTM_resCov,
                     nIter = 120000, burnIn = 20000, thin = 200, 
                     saveAt = paste0(mtmFolder, "/", treatment_i, '_'))
  )
  
  write.table(res_MTM$YHat, paste0(mtmFolder, "/", treatment_i, "_MTM_predictions.txt"), quote = F)
  write.table(res_MTM$resCov$R, paste0(mtmFolder, "/", treatment_i, "_MTM_resCov.txt"), quote = F)
  write.table(res_MTM$K[[1]]$G, paste0(mtmFolder, "/", treatment_i, "_MTM_G.txt"), quote = F)
  write.table(res_MTM$K[[1]]$U, paste0(mtmFolder, "/", treatment_i, "_MTM_U.txt"), quote = F)
  
  # diagnosis (log-likelihood)
  MTM_logLik <- read.table(paste0(mtmFolder, "/", treatment_i, '_logLik.dat'))
  png(paste0(mtmFolder, "/", treatment_i, "_logLik.png"))
  plot(MTM_logLik[[1]], cex = 0.5,
       main = paste('Log-likelihood, ', treatment_i),
       ylab = 'Log-likelihood')
  dev.off()
  
  # diagnosis (mu)
  MTM_mu <- readLines(paste0(mtmFolder, "/", treatment_i, '_mu.dat'))
  MTM_mu_1 <- sapply(MTM_mu, strsplit, split = ' ')
  names(MTM_mu_1) <- NULL
  MTM_mu_2 <- sapply(MTM_mu_1, function(x) {
    length_max <- max(sapply(MTM_mu_1, length))
    y <- rep(NA, length_max)
    y[1:length(x)] <- x
    as.numeric(y)
  })
  MTM_mu_3 <- t(as.data.frame(MTM_mu_2))
  png(paste0(mtmFolder, "/", treatment_i, "_mu.png"))
  matplot(MTM_mu_3, cex = 0.5, pch = 1:ncol(MTM_mu_3),
          main = paste0('mu, ', treatment_i, '_', the_year),
          ylab = 'mu')
  dev.off()
}
stopCluster(cl)

#### make corplot from "genetic covariance matrix" "residual covariance matrix" ######
#### make ggplot for the change of "obs_cor", "g_cor", "res_cor" #####
band_name <- c("dryWeight", "DTF")
condition <- c("W1", "W2", "W3", "W4")  


# make array to input the gCor between dryWeight and VIs
gCorDryWeightArray <- array(NA, dim = c(length(band_name), length(condition)))
dimnames(gCorDryWeightArray) <- list(band_name, condition)

# save the corplot as png
treatment <- c("C", "W5", "W10", "D")

# set working directory
data_folder <- mtmFolder

# read the cov files
g_cov_files <- list.files(data_folder, pattern = "MTM_G")
res_cov_files <- list.files(data_folder, pattern = "MTM_resCov")

g_cor_list <- list()
res_cor_list <- list()

png(paste0(mtmFolder, "/geneticCor.png"), 
    height = 1440, width = 1440, res = 216)
opar <- par(mfrow = c(2, 2))


for (file_index in 1:length(g_cov_files)) {
  # file_index <- 1
  treatmentEach <- treatment[file_index]
  g_cov_read <- g_cov_files[file_index]
  
  # make cor plot "g"
  g_cov_path <- paste0(data_folder, "/", g_cov_read)
  indexAll <- c("dryWeight", "DTF")
  cor_value_g0 <- MakeCorMat(g_cov_path, indexAll)
  cor_value_g <- cor_value_g0
  
  corrplot(cor_value_g, method = "shade", tl.srt = 45, addCoef.col = "black", 
           tl.col = "black", number.cex = .7)
  mtext(paste0("Genetic correlation ", treatmentEach)
        , at = 3, line = 2.5, cex = 1)
  
  gCorDryWeightArray[, file_index] <- cor_value_g["dryWeight", ]
} 
dev.off()
par(opar)
