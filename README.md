# Enamel_Sr

## Basic
Enamel_Sr support the data intake, processing, statistical analyses, and modeling components of the pre-print article entiled "Multi-substrate and micro-sampling analyses reveal spatial patterns of strontium isotope turnover in elephant molar enamel". The JAGS models made the following major improvements over the original [BITS framework](https://github.com/SPATIAL-Lab/BITS/tree/main): 1) it allows the model to be run more efficiently with rate scaling and bin thinning; 2) it allows the user to define the growth rate to sampling distance relationship to accommodate the logarithmic growth pattern of molar enamel and dentine.

## Data
The "data" folder contains all LA-ICP-MS, micromilled, and hand drilled enamel and dentine data used in this study, as well as laser log files and published data from [Yang et al., 2023](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14218) to support the analyses and plots associated with the article.

## Software requirements
The code was developed in R (ver. 4.3.1) calling the standalone JAGS (Just Another Gibbs Sampler) program, which is required before the code can be run. The JAGS program can be downloaded via [this link](https://sourceforge.net/projects/mcmc-jags/). Please make sure to download the version that is appropriate for your operating system.

To call JAGS from RStudio, the R packages "R2jags" are also required.

Other R packages specified in the file headers help to visualize data.

## Related publication and code
Yang et al., (2023). BITS: a Bayesian Isotope Turnover and Sampling model for strontium isotopes in proboscideans and its potential utility in movement ecology. Methods in Ecology and Evolution. https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14218

Yang, D. BITS: Bayesian Isotope Turnover and Sampling model. https://doi.org/10.5281/zenodo.7768584 

## How to cite this software
Yang, D. (2024), Enamel_Sr. https://github.com/Deming-Yang/Enamel_Sr