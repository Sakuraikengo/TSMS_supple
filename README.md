# TSMS_supple
Data which was used in the "Time-series Multi-spectral Imaging of Soybean for Predicting Biomass and Improving the Accuracy of Genomic Prediction" was uploaded here.

---

* scripts
    * 1.functionCode.R : register the function which was used in the scripts
    * 1.MTM_2.R : register the function of "MTM_2"
    * 2.0.makeAmat.R : make the additive relationship matrix
    * 2.1.0.spectralValueSel.R : calculate the vegetation index
    * 2.1.1.dataFolder.R : make the data.frame for the anarysis
    * 3.0.0.MTMFolder.R : make the multi-trait model(MTM)
    * 3.0.1.MTMFolderWithFlowerDate.R : add the flowering or not as fixed effect in "3.0.0"
    * 3.0.2.MTMFolderOnlyNDVI.R : make the multi-trait model using above-ground biomass (AGB) and only NDVI
    * 3.1.1.MTMpredictionWithFlowerDate.R : predict the AGB using the MTM
    * 3.1.2.MTMpredictionNDVI.R : predict the AGB using the MTM with only NDVI
    * 3.2.0.genomicPrediction.R : predict AGB using genomic information
    * 3.2.1.kernelPredictionWithFlowerDate.R : predict AGB using single/multi kernel model
    * 3.3.0.predictAllPooled.R : predict AGB over all treatments
    * 9.0.0.ANOVAdryweight.R : do the anova test for AGB data
    * 9.0.1.MTMFolderOnlyDTF.R : calculate the genetic correlation between days to flowering (DTF) and AGB
* 2019_M100_Xacti4eyeCamera_Images
    * field_data : This folder contains field data (ex. plot position, soil moisture).
    * MSdata : This folder contains multi-spectral data.
    * plantArea : This folder contains the plant area
    * raw : This folder contains spectral info each day (file size was too large to upload)
    * result : This folder contains the results of the each analysis
* figure
    Figures which were used in the paper
* genome
    Genome data (the file size of "Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz" was too large to upload)
* supple
    Supplementary data in the paper

