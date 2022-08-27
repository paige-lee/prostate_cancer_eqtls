## Loading libraries
library(data.table)
library(bigstatsr)

## Loading in the unzipped data
bedfile <- "/u/scratch/p/paigelee/Prostate_Cancer/genotypes/22.Prostate_Cancer.bed"
tempfile <- tempfile()
snp <- bigsnpr::snp_attach(bigsnpr::snp_readBed2(bedfile, backingfile = tempfile))

## Take transpose of snp$genotypes
snp$genotypes_transpose <- big_transpose(big_copy(snp$genotypes), backingfile = tempfile())

## Convert snp$genotypes_tranpose from an S4 object to a matrix
snp$genotypes_transpose <- snp$genotypes_transpose[]

## Add column names (sample ID) to the matrix snp$genotypes_tranpose
colnames(snp$genotypes_transpose) <- snp$fam$sample.ID

## Convert snp$genotypes_transpose from a matrix to a data.table
dim(snp$genotypes_transpose)
snp_data.table <- data.table("id" = snp$map$marker.ID, snp$genotypes_transpose)
head(snp_data.table)
