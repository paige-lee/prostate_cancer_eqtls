# prostate_cancer_eqtls

## Goal: 

Optimize a QTLtools pipeline for cis-eQTLs and trans-eQTLs to identify patterns of gene expression in prostate cancer patients


### 1. File Formatting

Convert input data files into the correct formats to use for MatrixEQTL in R. Formatting: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/html/Matrix_eQTL_main.html

### 2. Subsetting

Subsetted individuals in covariates, SNP location, and expression files so that all individuals are the same 

### 3. HCPs

HCPs = hidden covariate with prior knowledge

Wanted to obtain the HCPs from the expression matrix using Rhcpp (https://github.com/mvaniterson/Rhcpp)
