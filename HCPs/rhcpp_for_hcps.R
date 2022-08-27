# Load libraries
library(Rhcpp) 

######################################
# Subset the 500 covariates patients to match the (now) 348 expression patients
######################################
expression <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/mRNA/new_formatted_tumor_normal_Prostate_Cancer_Expression.csv", row.names = 1) # Read in expression data
tcga <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/Clinical/tcga.csv", row.names = 1) # Read in TCGA barcodes from covariates data

tcga <- as.character(unlist(tcga)) # Unlist
tcga <- toupper(tcga) # Convert to all uppercase
tcga <- gsub("-", ".", tcga) # Replace dashes with dots
length(tcga) # 490 patients in covariates data

expression_ids <- colnames(expression) # Get the TCGA barcodes of the expression data
expression_ids <- expression_ids[-1] # Remove the first column (ID) for now to avoid issues
expression_ids <- substr(expression_ids, 1, 12) # Remove unnecessary sample info so we can match up the TCGA barcodes
length(expression_ids) # 347 (without IDs column)

length(which(tcga %in% expression_ids)) # 347 samples from tcga that are in expression_ids
tcga <- tcga[which(tcga %in% expression_ids)] # Subset tcga IDs to just those in the expression data

covariates <- read.delim("/u/scratch/p/paigelee/Prostate_Cancer/Clinical/PRAD.clin.merged.txt") # Read in original covariates data
covariates_subset <- rbind(race = covariates[233, ], age = covariates[10, ], tcga = covariates[11, ]) # Bind the race, age, and TCGA covariates by row

white_or_NA_columns <- which(covariates_subset[1, ] == "white" | is.na(covariates_subset[1, ])) # Subset to get only the white patients or the NAs
length(white_or_NA_columns) # length is 490, which is 147 whites + 343 NAs
covariates_subset <- covariates_subset[, white_or_NA_columns]

covariates_ids <- covariates_subset[3, ] # TCGA barcodes for covariates data
covariates_ids <- as.character(covariates_ids) # Convert to a character vector
covariates_ids <- toupper(covariates_ids) # Convert to all uppercase
covariates_ids <- gsub("-", ".", covariates_ids) # Replace dashes with dots

length(which(covariates_ids %in% expression_ids)) # 347 of the patients from covariates are in the patients from expression
covariates_subset <- covariates_subset[, which(covariates_ids %in% expression_ids)]
dim(covariates_subset) # 3 covariates by 347 patients

covariates_subset <- covariates_subset[-1, ] # Remove race
covariates_subset <- covariates_subset[-2, ] # Remove TCGA

rownames(covariates_subset) <- NULL # Remove row names since we're not supposed to have them
covariates_subset <- cbind("id" = "age", covariates_subset) # Bind the "id" column to the covariates data

write.csv(covariates_subset, "/u/scratch/p/paigelee/Prostate_Cancer/Clinical/new_formattedCovariates.csv") # Write data
covariates <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/Clinical/new_formattedCovariates.csv", row.names = 1) # Read in data
dim(covariates) # 1 covariate by 348 patients (1 extra column)

######################################
# Use Nolan's normalization and standardization function on covariates data
######################################
geneNormalize <- function(vector) { # Function definition
  transformed <- log2(vector + 1) # Normalization
  (transformed - mean(transformed)) / sd(transformed) # Standardization
}

covariates <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/Clinical/new_formattedCovariates.csv", row.names = 1) # Read in covariates data
test <- as.numeric(covariates) # Convert covariates to an unnamed, numeric vector
test <- test[-1] # Remove the "age" entry from the covariates column
test <- geneNormalize(test) # Normalize and standardize the numeric entries 
rhcpp_covariates <- as.matrix(test) # Convert vector to a matrix
dim(rhcpp_covariates) # 347 by 1, which  matches the n (number of patients) by d (number of known covariates) requirement
write.csv(rhcpp_covariates, "rhcpp_covariates.csv") # Write this file
rhcpp_covariates <- read.csv("rhcpp_covariates.csv", row.names = 1) # Test it out
dim(rhcpp_covariates) # 347 by 1

######################################
# Use Nolan's normalization and standardization function on expression data
######################################
expression <- read.csv("/u/scratch/p/paigelee/Prostate_Cancer/mRNA/new_formatted_tumor_normal_Prostate_Cancer_Expression.csv", row.names = 1) # Read in expression data
test <- expression[, -1] # Remove non-numeric genes column
test <- t(test) # Take transpose of expression data since we need an n (number of patients) by g (number of genes) matrix
test <- apply(test, 2, geneNormalize) # Apply normalization and standardization to each column of the expression data 
length(which(is.na(test))) # 116,939 NAs
# 116,939 genes / 347 patients = 337 genes to remove

length(which(((rowSums(expression[, -1]) < 0.1) & (rowSums(expression[, -1]) >= 0)) | ((rowSums(expression[, -1]) < 347.1) & (rowSums(expression[, -1]) > 347.1))))
# I played with the above cutoff values so that there are 337 genes to remove
test <- expression[, -1]
test <- test[-which(((rowSums(expression[, -1]) < 0.1) & (rowSums(expression[, -1]) >= 0)) | ((rowSums(expression[, -1]) < 347.1) & (rowSums(expression[, -1]) > 347.1))), ]
dim(test) # 20165 genes by 347 patients
# 20502 genes - 20165 genes = 337 genes removed
test <- t(test) # Take transpose of expression data since we need an n (number of patients) by g (number of genes) matrix

test <- apply(test, 2, geneNormalize) # Apply normalization and standardization to each column of the expression data 
dim(test) # 347 patients by 20165 genes

write.csv(test, "rhcpp_expression.csv") # Write this file
rhcpp_expression <- read.csv("rhcpp_expression.csv", row.names = 1) # Test it out
dim(rhcpp_expression) # 347 by 20165

######################################
# Define r_hcp() function
######################################
r_hcp <- function(Z, Y, k, lambda1, lambda2, lambda3, iter=100) {
    ## convergence criteria
    tol <- 1e-6
    
    A <- matrix(0, ncol(Z), k)
    W <- matrix(0, nrow(Z), k)
    diag(W) <- 1
    n <- k*ncol(Y)
    ##B <- matrix(runif(n), ncol(Z), ncol(Y))
    B <- matrix((1:n)/n, k, ncol(Y))
    
    n1 <- nrow(Z)
    d1 <- ncol(Z)
    
    n2 <- nrow(Y)
    d2 <- ncol(Y)
    
    if(n1 != n2)
        message('number of rows in F and Y must agree')
    
    if (k < 1 | lambda1 < 0 | lambda2 < 0 | lambda3 < 0 )
        message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')
    
    ##predefine for slight preformance improvement
    diagB <- diag(k)
    diagW <- diag(k)
    diagA <- diag(nrow(A))
    U1 <- lambda1*solve(lambda1*crossprod(Z) + diagA*lambda3)%*%t(Z)
    
    if(iter > 0) {
        o <- numeric(iter)
        for(ii in 1:iter) {
            ##o[ii] <- norm(Y-W%*%B, type="F")^2 + lambda1*norm(W-Z%*%A, type="F")^2 + lambda2*norm(B, type="F")^2 + lambda3*norm(A, type="F")^2
            o[ii] <- sum((Y-W%*%B)^2) + sum((W-Z%*%A)^2)*lambda1 + sum(B^2)*lambda2 + lambda3*sum(A^2)           
            W <- (tcrossprod(Y, B) + lambda1*Z%*%A) %*% solve(tcrossprod(B) + lambda1*diagB)
            B <- solve(crossprod(W) + lambda2*diagW, crossprod(W,Y))            
            ##A <- solve(lambda1*crossprod(Z) + diagA*lambda3), lambda1*crossprod(Z,W))
            A <- U1%*%W
            
            if(ii > 1)  {
                if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
                    break
            }
            
        }
    }
    
    list(W=W, B=B, A=A, o=o, iter=ii)
}

######################################
# Use Rhcpp to get HCPs
######################################
covariates <- read.csv("rhcpp_covariates.csv", row.names = 1)
expression <- read.csv("rhcpp_expression.csv", row.names = 1)
covariates <- as.matrix(covariates)
expression <- as.matrix(expression)

r_hcp(Z = covariates , Y = expression, k = 5, lambda1 = 5, lambda2 = 5, lambda3 = 5)$W
