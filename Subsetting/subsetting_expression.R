# Load libraries
library(data.table)

# Read in TCGA barcodes from subsetted covariates file
tcga <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/Clinical/tcga.csv", row.names = 1)

# Reformat TCGA barcodes from covariates file to match the TCGA format of the expression data
tcga <- as.character(unlist(tcga)) # unlist
tcga <- toupper(tcga) # convert to all uppercase
tcga <- gsub("-", ".", tcga) # replace dashes with dots
length(tcga) # length is 490

# Read in formatted expression file
expression <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/mRNA/formatted_tumor_normal_Prostate_Cancer_Expression.csv")

# Reformat TCGA barcodes from expression data to match the TCGA format of the expression data
expression_ids <- substr(colnames(expression), 1, 12) # cut off the extra TCGA info about samples
expression_ids <- expression_ids[-1] # get rid of IDs column
length(expression_ids) # 498 patients in total

# How many TCGA barcodes match between covariates and expression data?
intersection <- RVenn::overlap(RVenn::Venn(list(expression_ids, tcga))) # intersection between expression IDs and covariates
indices <- which(expression_ids %in% intersection) # Which entries from expression_ids are in the intersection?
length(indices) # 489 entries from expression_ids are in intersection

# Subset expression data based on intersection with covariates
test <- expression[, -1] # Remove ID column for now since it's the first column, but the intersection doesn't include the ID column, so the index 1 will be different 
test <- test[, indices] # Only the columns that match the indices
test <- data.table("id" = expression[, 1], test) # Add back the first column ID

# Delete the duplicated TCGA  barcode (duplicated patient)
expression_ids <- colnames(test) # Expression IDs are the column names of the expression matrix
expression_ids <- expression_ids[-1] # Remove the first column (ID)
expression_ids <- substr(expression_ids, 1, 12) # Remove uneccessary sample information
length(which(duplicated(expression_ids) == TRUE)) # There is 1 duplicate
which(duplicated(expression_ids) == TRUE) # 372 is the duplicate
which(expression_ids == expression_ids[372]) # Which entry equals the duplicate entry? 371 and 372
expression_ids[371] # "TCGA.V1.A9O5"
expression_ids[372] # "TCGA.V1.A9O5"
expression_ids <- colnames(test) # Expression IDs are the column names of the expression matrix
expression_ids <- expression_ids[-1] # Remove the first column (ID)
expression_ids[371] # Without truncating the TCGA barcode, this is TCGA.V1.A9O5.06A.11R.A41O.07
expression_ids[372] # Without truncating the TCGA barcode, this is TCGA.V1.A9O5.01A.11R.A41O.07
# The duplicates are really different samples from the same patient
# Arjun said to keep the 06A sample (entry 371) but delete the 01A sample (entry 372)
colnames(test)[373] # TCGA.V1.A9O5.01A.11R.A41O.07 confirmed (372 + 1 for the ID column)
test[, 373] # TCGA.V1.A9O5.01A.11R.A41O.07 confirmed
test <- test[, -373] # Delete TCGA.V1.A9O5.01A.11R.A41O.07
dim(test) # 20502 by 489 (1 extra column for IDs, so 488 patients)

# Write the new formatted file and test it out
write.csv(test, "new_formatted_tumor_normal_Prostate_Cancer_Expression.csv")
test <- read.csv("new_formatted_tumor_normal_Prostate_Cancer_Expression.csv", row.names = 1)
dim(test) # 20502 by 489 (one extra column for IDs)

# Subset expression file (488 patients) for the patients in the genotype files (348 patients)
expression <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/mRNA/new_formatted_tumor_normal_Prostate_Cancer_Expression.csv", row.names = 1) # Read in expression data
genotype <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/genotypes/21.Prostate_Cancer_genotype.csv", header = TRUE) # Read in genotype data

dim(expression) # 20502 genes by 489 samples (one extra column for IDs)
dim(genotype) # 58823 SNPs by 348 samples (one extra column for IDs)

expression_cols <- colnames(expression) # Column names of expression data are the samples
genotype_cols <- colnames(genotype) # Column names of genotype data are the samples

expression_cols <- substr(expression_cols, 1, 12) # Remove extra information so the IDs match
genotype_cols <- substr(genotype_cols, 1, 12) # Remove extra informationn so the IDs match

length(which(expression_cols %in% genotype_cols)) # 348 entries of expression_cols are in genotype_cols

expression <- expression[, which(expression_cols %in% genotype_cols)] # Subset expression data for only samples genotype data
dim(expression) # 20502 genes by 348 samples (one extra column for IDs)

expression <- data.table(expression) # Convert to a data.table

# Write the new formatted file and test it out
write.csv(expression, "/u/scratch/p/paigelee/Prostate_Cancer/mRNA/new_formatted_tumor_normal_Prostate_Cancer_Expression.csv") # Write data
expression <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/mRNA/new_formatted_tumor_normal_Prostate_Cancer_Expression.csv", row.names = 1) # Read data
dim(expression) # 20502 by 348 (one extra column for IDs)
