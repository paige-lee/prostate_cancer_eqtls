# Load libraries
library(data.table)
library(genio)

# Load data
bimfile <- read_bim("22.Prostate_Cancer.bim")
head(bimfile)

# Construct formatted data.table
data.table <- data.table(snp = bimfile$id, chr = bimfile$chr, pos = bimfile$pos)
head(data.table)
