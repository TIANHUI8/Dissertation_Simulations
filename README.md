# Dissertation_Simulations
A set of R scripts used to simulate GWAS data in order to test variable selection methods

## Contents
### bird analysis.txt

This is the R code for variable selection analysis for bird data in a txt file with 6 methods:
* single marker regression
* lasso
* ridge regression
* elastic-net
* Bayesian ridge regression
* BayesA
This script also include the scatterplot for orderings of the 7 methods
### BAP_analysis.R
This script will analysis the sorghum data with 7 methods above in bird data;
Since the number of SNPs for the sorghum data is much larger than bird data, this script will run on Cluster

### Sim_a.R
This script will simulate SNP marker data with 5 different possible confounding factors:
* n (sample size) vs. p (marker number) ratio
* correlation among markers
* family structure
* measurement error
* model complexity

This script can create 4 simulated data sets with different combinations of the above factors.
