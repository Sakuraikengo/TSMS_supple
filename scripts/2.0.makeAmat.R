library(gaston)
library(RAINBOWR)

# make the result folder
resultFolder <- "2019_M100_Xacti4eyeCamera_Images/result"
if (!file.exists(resultFolder)) {
  dir.create(resultFolder)
}

# set the seed
seedIndCsv <- paste0(resultFolder, "/seedInd.csv")
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(paste0(resultFolder, "/seedInd.csv"), row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:500, 10, replace = F)
  write.csv(x = seedInd, file = paste0(resultFolder, "/seedInd.csv"))
}


# read the vcf file
genomeFile <- "genome/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz"
genomeRaw <- read.vcf(genomeFile)

# filtering the SNPs data
genomeFil <- select.snps(genomeRaw, condition = maf >= 0.05)
genomeFil <- LD.thin(genomeFil, threshold = 0.8)

genomeMat <- as.matrix(genomeFil) - 1
amat <- calcGRM(genoMat = genomeMat, methodGRM = "A.mat")
write.csv(x = amat, file = "genome/amat173583SNP.csv")
