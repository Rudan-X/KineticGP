# KineticGP
KineticGP is an computational framework combining [kinetic model of C4 photosynthesis](https://github.com/yuwangcn/C4_dynamic_model) with genomic prediction.
It allows the identification of  key kinetic parameters identified via sensitivity analysis 

## OS

- Code was tested on Debian 6.1.128-1 and Windows 10



## Dependencies

- Matlab (tested with version 2023b)
- R (tested with version 4.3.1)
- [PESTO](https://github.com/ICB-DCM/PESTO/tree/master)



## Sensitivity analysis of kinetic parameters
To estimate kinetic parameters of top control coefficient for each genotype, run the code [optim_MCMC_sampling.m](https://github.com/Rudan-X/KineticGP/blob/master/code/matlab/parameterization/optim_MCMC_sampling.m) 

To perform parallel estimation of all genotypes on linux HPC, submit batch [matlab_sensitivity.job](https://github.com/Rudan-X/KineticGP/blob/master/code/bash/matlab_parameterization.job) to Slurm 


## Estimation of kinetic parameters
Please run the code KineticGP\code\matlab\parameterization\parallel_parameterization.m or submit KineticGP\code\bash\matlab_parameterization.job on linux HPC


## Genomic Prediction of kinetic parameters and photosynthetic rate
Please run the codes in KineticGP\code\R\GenomicPrediction 


## Reproduce the results presented in the paper

Please run the codes in folder KineticGP\code\matlab\generateFigures and KineticGP\code\R\generateFigures

