# Load libraries
library(data.table)

# Read in covariates data
covariates <- read.delim("PRAD.clin.merged.txt")
covariates[, 1] # View all covariate names to see which ones I want
dim(covariates) # 1459 covariates by 500 columns 

# Bind the race, age, gender, and tumor level covariates by row
covariates_subset <- rbind(race = covariates[233, ], age = covariates[10, ], gender = covariates[205, ], tumor_level = covariates[334, ])
dim(covariates_subset) # 4 covariates by 500 columnns

# Break down of patients by race
table(as.character(covariates_subset[1, ]))
length(which(is.na(covariates_subset[1, ])))
# asian: 2, black or african american: 7, patient.race_list.race: 1 (this is the covariates column, so ignore), white: 147, NAs: 343
# Sum: 2 + 7 + 147 + 343 = 499 patients total, which matches 500 columns - 1 covariates column

# Subset to get only the white patients or the NAs
white_or_NA_columns <- which(covariates_subset[1, ] == "white" | is.na(covariates_subset[1, ]))
length(white_or_NA_columns) # length is 490, which is 147 whites + 343 NAs
covariates_subset <- covariates_subset[, white_or_NA_columns]
dim(covariates_subset) # 4 covariates by 490 columnns

# Convert from data.frame to data.table
covariates_subset <- data.table(covariates_subset)

# Add first column with covariate names
covariates_subset <- data.table("id" = c("race", "age", "gender", "tumor_level"), covariates_subset)

# Check which rows/covariates have NAs
any(is.na(covariates_subset[1, ])) # TRUE (there are NAs for race)
any(is.na(covariates_subset[2, ])) # FALSE (no NAs for age)
any(is.na(covariates_subset[3, ])) # FALSE (no NAs for gender)
any(is.na(covariates_subset[4, ])) # TRUE (there are NAs for tumor level)

# How many NAs does tumor level have?
length(which(is.na(covariates_subset[4, ]))) # 140 out of 490 are NAs, which is too many NAs for a covariate

# Remove race and gender since all patients are white or NA and male, and remove tumor_level since there are too many NAs
covariates_subset <- covariates_subset[-1, ] # Remove race
covariates_subset <- covariates_subset[-2, ] # Remove gender
covariates_subset <- covariates_subset[-2, ] # Remove tumor level
dim(covariates_subset) # 1 covariate by 491 columns (first column contains covariate names)

# Write the file
write.csv(covariates_subset, "formattedCovariates.csv")
