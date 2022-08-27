# Load the libraries
library(data.table)

# Read in the data
test <- read.delim("tumor_normal_Prostate_Cancer_Expression.txt")
test <- test[, 1]
test <- as.character(test)

# Read in the annotation file (subset with relevant columns)
subset <- read.csv("subset.csv")
subset <- subset[, c(2, 3, 4, 5)]

# Extract the genes, matches, and subset the annotation file based on the data
subset_genes <- as.character(subset$gene_name)
subset_matches <- which(subset_genes %in% test)
subset_genes2 <- subset_genes[subset_matches]
subset2 <- subset[subset_matches, ]

# Create the output data.table file with default values
output <- data.table("geneid" = test, "chr" = rep("", length(test)), "s1" = rep(0, length(test)), "s2" = rep(0, length(test)))

print("Ready to start matching")

# Going through each gene in from the data to find its match (if any) in the annotation file
for (i in 1:length(test)) {
  for (j in 1:length(subset_genes2)) {
    if (identical(test[i], subset_genes2[j])){ # If the gene names are identical...
      output[i, 2] <- as.character(subset2[j, 1]) # chr in the new data.table becomes chr from the annotation file 
      output[i, 3] <- as.numeric(subset2[j, 2]) # s1 in the new data.table becomes start from the annotation file
      output[i, 4] <- as.numeric(subset2[j, 3]) # s2 in the new data.table becomes end from the annotation file
    }
  }
}

warnings()

print("Finished matching")

# Write the new data.table as a .csv file
write.csv(output, "tumor_normal_Prostate_Cancer_gene_location.csv")
