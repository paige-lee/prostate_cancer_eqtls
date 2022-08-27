# Import libraries
library(biomaRt)

# Read in data
data <- read.csv("normal_Prostate_Cancer_gene_location.csv")
data <- data[, 2:5]
missing <- which(data$chr == "")

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

print("Ready to start matching")

for (i in 1:length(missing)) {
temp <- getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position", "end_position"), filters = "hgnc_symbol", values = as.character(data[missing[i], 1]), mart = human)

warnings()

if (dim(temp) == c(0, 4)) { # If the missing gene is not found in this dataset, then leave the info blank as is and move on
next
}

else { # If the missinng gene is found in this dataset, then fill in the missing info
data[missing[i], 2] <- as.character(temp$chromosome_name) # chromosome 
data[missing[i], 3] <- as.numeric(temp$start_position) # start position
data[missing[i], 4] <- as.numeric(temp$end_position) # end position
}

warnings()

}

print("Ready to start writing")

write.csv(data, "normal_Prostate_Cancer_gene_location2.csv")
