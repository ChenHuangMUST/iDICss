**iDICss robustly predicts melanoma immunotherapy response by synergizing genomic and transcriptomic knowledge via independent component analysis**

Repository for the paper "Integration of genomic and transcriptomic data using independent component analysis to predict immunotherapy response in melanoma"

This repository contains all the scripts for the analyses presented in the paper. This includes scripts for running ICA, obtaining immune-related ICs, calculating the iDIC matrix and constructing the iDIC scoring system.

**Contents**

 - runICA/: This folder contains the scripts for running ICA.
 - getimmIC.R: Filtering for independent components related to immunology.
 - getiDICmatrix.R: Calculation of TIME-driver independent components profile by integrating signature components and somatic mutation data.
 - getiDICss.R: Development of iDIC-based risk score system.
 - getMLcombine.R：A function that develops robust predictive models with excellent performance, using patient prognosis as the outcome.
 - plot.heatmap.R：Plot the heatmap of the C-index of all combined prediction models.
 - compare.performance.R：Calculating the performance measurements of all combined prediction models.

**Installation**

All packages can be installed via Bioconductor (https://www.bioconductor.org/), GitHub, and CRAN. Generally, a couple of minutes is needed for installing each package.

**Contact**

For any questions or feedback, please reach out to Chen Huang at chuang@must.edu.mo.


