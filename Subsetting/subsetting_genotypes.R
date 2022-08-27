# Mostly Nolan's code except this is the intersection between covariates and clinical since I already subsetted the expression file
# Load libraries
library(bigstatsr)
library(RVenn)

# Load data (covariates data, which is subsetted for whites and NAs, and one chromosome of genotype data)
genotype <- bigsnpr::snp_attach(bigsnpr::snp_readBed2("22.Prostate_Cancer.bed"))
clinical <- data.table::fread("/u/scratch/p/paigelee/Prostate_Cancer/Clinical/tcga.csv", header = TRUE)

# Manipulate TCGA bar codes so that their formats are the same
p_genotype <- toupper(genotype$fam$sample.ID)
p_genotype_sub <- substr(p_genotype, 1, 12)

p_clinical <- toupper(as.character(unlist(clinical[, 2])))
p_clinical_sub <- substr(p_clinical, 1, 12)

# Find the intersection of TCGA bar codes between covariates and genotype data (will select just whites and NAs in genotype data)
intersection <- RVenn::overlap(RVenn::Venn(list(p_clinical_sub, p_genotype_sub)))

# Write the TCGA bar codes as a .csv file to use in genotype file formatting
data.table::fwrite(data.frame(intersection), "uniquePatients.csv")

# Remove unnecessary files
file.remove(list.files()[stringr::str_detect(list.files(), ".rds")])
file.remove(list.files()[stringr::str_detect(list.files(), ".bk")])
