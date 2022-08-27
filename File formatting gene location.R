# Load libraries
library(data.table)

# Load data
txtfile <- read.delim("Prostate_Cancer_CNV.txt")
head(txtfile)

# Construct formatted data.table
data.table <- data.table("Sample" = txtfile$Sample, "chr" = txtfile$Chromosome, "s1" = txtfile$Start, "s2" = txtfile$End)

# Issue: not sure how to get gene ID from sample for the 1st column