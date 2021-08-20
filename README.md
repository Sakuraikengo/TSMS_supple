# TSMS_supple
Data which was used in the "Time-series Multi-spectral Imaging of Soybean for Predicting Biomass and Improving the Accuracy of Genomic Prediction" was uploaded here.

If you want to watch the "supplementary_data", please click the "supplementary_data.zip". Then, Download the zip file. 

---

* scripts
    * 1.functionCode.R : register the function which was used in the scripts
    * 1.MTM_2.R : register the function of "MTM_2"
    * 2.0.makeAmat.R : make the additive relationship matrix
    * 2.1.0.spectralValueSel.R : calculate the vegetation index
    * 2.1.1.dataFolder.R : make the data.frame for the anarysis
    * 3.0.0.MTMFolder.R : make the multi-trait model(MTM)
    * 3.0.1.MTMFolderWithFlowerDate.R : add the flowering or not as fixed effect in "3.0.0"
    * 3.1.0.MTMprediction.R : predict the above-ground biomass (AGB) using the MTM
    * 3.1.1.MTMpredictionWithFlowerDate.R : add the flowering or not as fixed effect in "3.1.0"
    * 3.2.0.kernelPrediction.R : predict AGB using single/multi kernel model
    * 3.2.1.kernelPredictionWithFlowerDate.R : add the flowering or not as fixed effect in "3.2.0"
    * 3.3.0.predictAllPooled.R : predict AGB over all treatments
* supplemenatry_data.zip
    * 2019_M100_Xacti4eyeCamera_Images
        * field_data : This folder contains field data (ex. plot position, soil moisture).
        * MSdata : This folder contains multi-spectral data.
        * plantArea : This folder contains the plant area
        * result : This folder contains the results of the each analysis
    * figure
    Figures which were used in the paper
    * genome
    Genome data 
    * supple
    Supplementary data in the paper

