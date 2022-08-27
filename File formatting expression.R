# Load libraries
library(data.table)

# Load data
txtfile <- read.delim("normal_Prostate_Cancer_Expression.txt")
head(txtfile)

# Construct formatted data.table
dim <- dim(txtfile)
data.table <- data.table("id" = txtfile$Gene, txtfile[,2:dim[2]])
head(data.table)