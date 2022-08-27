#################################
# Issue 1: Nolan's gene location file has ensembl IDs, but I need hgnc symbols for my data
#################################

# Load libraries
library(biomaRt)
library(RVenn)

# Use biomaRt to find the corresponding hgnc symbols to the ensembl IDs in the gene location file
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # Load the hsapiens_gene_ensembl dataset

df <- read.csv("formatted_gene_loc.csv") # Read in gene location data
df <- df$id # Extract only the ensembl IDs
df <- as.character(df) # Convert to a character vector

G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = df, mart = mart) 
# Use the ensembl IDs to get the corresponding hgnc symbols 

write.csv(G_list, "ensembl_to_hgnc.csv") # Write this conversion as a .csv file

# Replace hgnc symbols with ensembl IDs in the gene location file
conversion <- read.csv("ensembl_to_hgnc.csv", row.names = 1) # Read in conversion data
dim(conversion) # Is 21261 genes by 2 columns (ensembl, hgnc)
gene_location <- read.csv("formatted_gene_loc.csv") # Read in gene location data
dim(gene_location) # Is 23651 genes by 4 columns (ID, chr, s1, s2)
# There's a difference of 23651 (gene location) - 21261 (conversion) = 2390 genes that don't have an hgnc symbol

gene_location[, 1] <- ifelse(gene_location[, 1] %in% conversion[, 1], conversion[, 2], NA)
# Replace the ensembl IDs in the gene location file with the hgnc symbols in the conversion file if the ensembl IDs match
# or with an NA otherwise

write.csv(gene_location, "hgnc_gene_loc.csv") # Write the gene location file
test <- read.csv("hgnc_gene_loc.csv", row.names = 1) # Test it out

#################################
# Issue 2: not all of the ensembl IDs were able to be converted to hgnc symbols, and not all of the genes in the expression data were found in the 
# reformatted gene location file
# Solution: use the intersection of genes from the gene location file and the expression data (delete genes that don't intersect)
#################################

# Read in data
gene_loc <- read.csv("hgnc_gene_loc.csv", row.names = 1) 
expression <- read.csv("new_formatted_tumor_normal_Prostate_Cancer_Expression.csv", row.names = 1) 

# Remove NAs (genes that couldn't be converted from ensembl to hgnc) from gene location file 
remove <- which(is.na(gene_loc[, 1])) # Indices of NA genes to remove
gene_loc <- gene_loc[-remove, ] # Update gene location to not include NA genes
dim(gene_loc) # 21615 by 4

# Find the intersection of genes between gene location and expression
intersection <- RVenn::overlap(RVenn::Venn(list(gene_loc[, 1], expression[, 1])))
length(intersection) # 15445 genes in intersection 
# Gene location: 21615 - 15445 = 6170 genes in gene location that are not in expression 
# Expression: 20502 - 15445 = 5057 genes in expression that are not in gene location

# Which genes should we keep in each file?
indices_gene_loc <- which(gene_loc[, 1] %in% intersection) # Which genes from gene location should we keep (in intersection)?
indices_expression <- which(expression[, 1] %in% intersection) # Which genes from expression should we keep (in intersection)?

length(indices_gene_loc) # 17927 we don't want this (larger than intersection)
length(indices_expression) # 15445 this is what we want (size of intersection)
length(which(duplicated(gene_loc[indices_gene_loc, 1]) == TRUE)) # 2482 duplicates, and 17927 - 2482 = 15445 so perfect
# We can keep the duplicates for now and have 17927 genes in gene location and 15445 genes in expression

gene_loc <- gene_loc[indices_gene_loc, ] # Only keep genes in the intersection
expression <- expression[indices_expression, ] # Only keep genes in the intersection
dim(gene_loc) # 17297 by 4
dim(expression) # 15445 by 348

# Write the new data
write.csv(gene_loc, "hgnc_gene_loc.csv")
write.csv(expression, "new_formatted_tumor_normal_Prostate_Cancer_Expression.csv")

#################################
# Issue 3: The gene names in gene location and expression need to be exactly the same and no duplicates
#################################

# Investigate the duplicates issue
gene_loc <- read.csv("hgnc_gene_loc.csv", row.names = 1) # Read in data
duplicates <- which(duplicated(gene_loc[, 1]) == TRUE) # Which rows/genes are duplicates?
length(duplicates) # There are 2482 duplicated genes
gene_loc[duplicates, ] # What do the duplicated rows look like? The gene names are duplicated, but their chr, s1, and s2 information are different
# This discrepency probably happened when I converted from ensembl to hgnc
# Goal: add chr, s1, and s2 information to the query to figure out what the correct information is for each gene and then delete the rows of duplicated genes

# Use biomaRt to find out what's the correct information
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # Load the hsapiens_gene_ensembl dataset
listFilters(mart) # What are the names of the filters we want to use?
listAttributes(mart) # What are the names of the attributes we want to use?

df <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/mRNA/formatted_gene_loc.csv") # Read in Nolan's gene location file
df <- df$id # Extract just the ensembl name of the gene
df <- as.character(df) # Convert to a character vector

G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), values = df, mart = mart)
# Filter by ensembl ID to obtain the ensembl ID, hgnc symbol, chromosome name, start position, and end position
dim(G_list) # 20733 by 5

#################################
# Issue 4: We weren't using the most recent version of Nolan's gene location file
#################################

# Use biomaRt to find the corresponding hgnc symbols to the ensembl IDs in the gene location file
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) # Load the hsapiens_gene_ensembl dataset

df <- read.csv("formatted_gene_loc.csv") # Read in gene location data
df <- df$id # Extract only the ensembl IDs
df <- as.character(df) # Convert to a character vector

G_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = df, mart = mart) 
# Use the ensembl IDs to get the corresponding hgnc symbols 

write.csv(G_list, "ensembl_to_hgnc.csv") # Write this conversion as a .csv file

# Replace hgnc symbols with ensembl IDs in the gene location file
conversion <- read.csv("ensembl_to_hgnc.csv", row.names = 1) # Read in conversion data
dim(conversion) # Is 20733 genes by 2 columns (ensembl, hgnc)
gene_location <- read.csv("formatted_gene_loc.csv") # Read in gene location data
dim(gene_location) # Is 22909 genes by 4 columns (ID, chr, s1, s2)
# There's a difference of 22909 (gene location) - 20733 (conversion) = 2176 genes that don't have an hgnc symbol

gene_location[, 1] <- ifelse(as.character(gene_location[, 1]) %in% as.character(conversion[, 1]), as.character(conversion[, 2]), NA)
# Replace the ensembl IDs in the gene location file with the hgnc symbols in the conversion file if the ensembl IDs match
# or with an NA otherwise

write.csv(gene_location, "hgnc_gene_loc.csv") # Write the gene location file
test <- read.csv("hgnc_gene_loc.csv", row.names = 1) # Test it out
