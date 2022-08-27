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

### 4. Using matrix eQTL

My code for this part accidentally got deleted. I made a script that incorporates the files from the previous steps to pass through the matrix eQTL function. Then, I made a bash script that executes the R script and submit it as a job to the queue of the computing cluster.
