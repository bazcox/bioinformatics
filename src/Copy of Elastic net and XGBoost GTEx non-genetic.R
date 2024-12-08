# GTEx
# Elastic net 
# standard non-oversampled test set
# with non-genetic data

###################################
#                                 #
#      Probe 1 - cg04692870       #
#                                 #            
###################################
library(janitor) 
library(psych)
library(tidyverse)
library(xgboost)
library(DataExplorer)
library(caret)
library(randomForest)
library(readxl)
library(qdapRegex)
library(factoextra)
library(mgcv)
library(olsrr)
library(pls)
library(data.table)
library(leaps)
library(MASS)
library(neuralnet)
library(iml)
library(missForest)
library(Boruta)
library(glmnet)
library(ggplot2)
library(readxl)
library(dplyr)

library(readxl)

GenotypeFile <- read.csv("GTEx08snps_MAF10_genotype.csv")
pca <- read.csv("probe1_v2_genotype_subset_LDP_pca.eigenvec_5PCs_allsamples.csv")
demo <- read.csv("GUSTO_nongenetic_demo_withNA.csv")

Probe1 <- read.csv("Probe1_cg04692870_v2.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")


#merge the files by FID
MainData <- GenotypeFile
MainData <- merge(MainData, demo, by = "FID")
MainData <- merge(MainData, pca, by = "FID")

MainData <- merge(MainData, Probe1, by = "FID")


test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#replace 0 with NA
train.set <- na_if(train.set, 0)
test.set <- na_if(test.set, 0)

# mode imputation
# highest no. of N/A in a col: 51
my_mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  if (length(xmode) > 1) xmode <- sample(xmode, 1)
  return(xmode)
}

for (var in 1:ncol(train.set)) {
  if (class(train.set[,var]) %in% c("character", "factor")) {
    train.set[is.na(train.set[,var]),var] <- my_mode(train.set[,var], na.rm = TRUE)
  }
}

for (var in 1:ncol(test.set)) {
  if (class(test.set[,var]) %in% c("character", "factor")) {
    test.set[is.na(test.set[,var]),var] <- my_mode(test.set[,var], na.rm = TRUE)
  }
}

# mode imputation for non-genetic cols
calc_mode <- function(x){
  
  # List the distinct / unique values
  distinct_values <- unique(x)
  
  # Count the occurrence of each distinct value
  distinct_tabulate <- tabulate(match(x, distinct_values))
  
  # Return the value with the highest occurrence
  distinct_values[which.max(distinct_tabulate)]
}


train.set <- train.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))
test.set <- test.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))


#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

#mean imputation of mother age recruitment
train.set$mother_age_recruitment[is.na(train.set$mother_age_recruitment)] <- mean(train.set$mother_age_recruitment, na.rm = TRUE)
test.set$mother_age_recruitment[is.na(test.set$mother_age_recruitment)] <- mean(test.set$mother_age_recruitment, na.rm = TRUE)



### step 1
# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


### step 2
# for making sure train.set always has at least 2 samples for each category/ class for each column
# for columns that have a category/ class that has only 1 sample (eg X.83 only has 1 G G)
# add that category (eg G G) to that column in a new row
# and for the rest of the columns: add a random genotype already existing in the column

# Count the number of samples in each category class for each column
for (col in names(train.set)) {
  counts <- table(train.set[[col]])
  
  # Check if there is any category with only one sample
  if (sum(counts == 1) > 0) {
    # Print the column name if there is at least one category with one sample
    print(col)
  }
}

# Create a new row with all missing values
# ignore if no snps output from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 562))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 191] <- "A A" #190 + 1
new_row[, 192] <- "G G" #191 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(191,192)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:562, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(549), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



# Loop through the columns with 2 classes
# ignore if no output in the variable from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(test.set)))

for (col in cols_with_2_classes) {
  print(col)
  # Find the missing class in df1
  missing_class <- setdiff(levels(factor(train.set[[col]])), levels(factor(test.set[[col]])))
  
  
  #colnames(new_row) <- colnames(test.set)
  new_row[[col]] <- missing_class
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add the new row to test.set
# ignore if no output from prev code
colnames(new_row) <- colnames(test.set)
test.set <- rbind(test.set, new_row)

## for checking
# Make sure all columns in test.set have 3 unique values
for (col in names(test.set)) {
  if (length(unique(test.set[[col]])) < 3) {
    stop("Not all columns in test.set have 3 unique values.")
  }
}


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- seq (0, 1, by = 0.05) # 0.00 0.05 0.10 0.15 0.20 0.25 etc
lambda.vec <- 0:30 # 0  1  2  3  4  5  6 etc

# expand.grid: 
# Create a Data Frame from All Combinations of Factor Variables
elastic.grid <- expand.grid (alpha = alpha.vec, lambda = lambda.vec)
k = 10
repeats = 5

# repeatedcv:
# repeatedly performs X-fold cross-validation on the training data, 
# i.e. if you specify 5 repeats of 10-fold cross-validation, 
# it will perform 10-fold cross-validation on the training data 5 times, 
# using a different set of folds for each cross-validation
tr.Control <- trainControl (method = "repeatedcv",
                            number = k,
                            repeats = repeats,
                            search = "grid")


set.seed (1403)
elastic <- train (cg04692870 ~ ., data = train.set, 
                  method = 'glmnet', 
                  trControl = tr.Control, 
                  verbose = FALSE,
                  tuneGrid = elastic.grid)

elastic

elastic_trainresults <- data.frame(matrix(ncol = 4, nrow = 0))

elastic_trainresults <- cbind (elastic$results$alpha, elastic$results$lambda, elastic$results$RMSE, elastic$results$Rsquared)

colnames(elastic_trainresults) <- c('alpha', 'lambda', 'RMSE', 'Rsquared')

elastic_trainresults <- data.frame (elastic_trainresults)

# alpha = 0: ridge regression
# alpha = 1: lasso regression

plot (elastic)

lambda.opt <- elastic$bestTune$lambda
alpha.opt <- elastic$bestTune$alpha

index.opt <- which(elastic$results$lambda == lambda.opt  & 
                     elastic$results$alpha == alpha.opt )

results.elastic.opt <- elastic$results[index.opt, ]

results.elastic.opt


#predict on test.set without last row of dummy variables
elastic.pred <- predict(elastic, newdata = test.set[,1:561])


df <- data.frame(RMSE.el = RMSE(elastic.pred, test.set$cg04692870),
                 Rsquare.el = caret::R2(elastic.pred, test.set$cg04692870))



# XGBoost
#-----------------------------------------------------------------------------------------------
eval_results <- function (true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(RMSE = RMSE,
             Rsquare = R_square)
}
#-----------------------------------------------------------------------------------------------

# setting up xbg.DMatrix

# convert to sparse matrix
sparse_matrix <- sparse.model.matrix(cg04692870 ~ ., data = train.set)[,-1]
dtrain_mQTL = xgb.DMatrix(sparse_matrix, label = as.matrix(train.set$cg04692870))

### parameter tuning 100 rounds, random selection of parameters
best.param = list()
best.cv = Inf
best.cv.index = 0
best.nround = 0
set.seed(123)

for (iter in 1:50) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = sample(1:10, 1),
                eta = runif(1, 0.001, 0.05),
                subsample = runif(1, .5, .9),
                colsample_bytree = runif(1, .5, .9) 
  )
  cv.nround = 1000
  #uncomment nthread=16 if you want parallel processing, and you know the number of cores you have, i have 16 so i use 16 here
  mdcv <- xgb.cv(data=dtrain_mQTL, params = param, 
                 # nthread=6,
                 nfold=10, nrounds=cv.nround,
                 # verbose = FALSE, 
                 early_stopping_rounds=100)
  
  min.cv = min(mdcv$evaluation_log$test_rmse_mean)
  min.cv.index = which.min(mdcv$evaluation_log$test_rmse_mean)
  
  if (min.cv < best.cv) {
    best.cv = min.cv
    best.cv.index = min.cv.index
    best.param = param
  }
}

# best model 

best.param
nround = best.cv.index
nround

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)
xgb.save(best.boost, "xgboost.GTEx.model.non-genetic.probe1")


### predictions
xgboost_train_predictions <- as.data.frame(predict(best.boost, sparse_matrix))   

# convert test set to sparse matrix
sparse_matrix_test <- sparse.model.matrix(cg04692870 ~ ., data = test.set)[,-1]

xgboost_predictions <- as.data.frame(predict(best.boost, sparse_matrix_test))
test.set$cg04692870
xgboost_predictions

eval_results (train.set$cg04692870, xgboost_train_predictions, train.set) #train
eval_results (test.set$cg04692870, xgboost_predictions, test.set) #test


###################################
#                                 #
#      Probe 3 - cg09322432       #
#                                 #            
###################################
library(janitor) 
library(psych)
library(tidyverse)
library(xgboost)
library(DataExplorer)
library(caret)
library(randomForest)
library(readxl)
library(qdapRegex)
library(factoextra)
library(mgcv)
library(olsrr)
library(pls)
library(data.table)
library(leaps)
library(MASS)
library(neuralnet)
library(iml)
library(missForest)
library(Boruta)
library(glmnet)
library(ggplot2)
library(readxl)
library(dplyr)

library(readxl)

GenotypeFile <- read.csv("GTEx08snps_MAF10_genotype.csv")
pca <- read.csv("probe3_genotype_subset_LDP_pca.eigenvec_5PCs.csv")
demo <- read.csv("GUSTO_nongenetic_demo_withNA.csv")

Probe1 <- read.csv("Probe3_cg09322432.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")


#merge the files by FID
MainData <- GenotypeFile
MainData <- merge(MainData, demo, by = "FID")
MainData <- merge(MainData, pca, by = "FID")

MainData <- merge(MainData, Probe1, by = "FID")


test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#replace 0 with NA
train.set <- na_if(train.set, 0)
test.set <- na_if(test.set, 0)

# mode imputation
# highest no. of N/A in a col: 51
my_mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  if (length(xmode) > 1) xmode <- sample(xmode, 1)
  return(xmode)
}

for (var in 1:ncol(train.set)) {
  if (class(train.set[,var]) %in% c("character", "factor")) {
    train.set[is.na(train.set[,var]),var] <- my_mode(train.set[,var], na.rm = TRUE)
  }
}

for (var in 1:ncol(test.set)) {
  if (class(test.set[,var]) %in% c("character", "factor")) {
    test.set[is.na(test.set[,var]),var] <- my_mode(test.set[,var], na.rm = TRUE)
  }
}

# mode imputation for non-genetic cols
calc_mode <- function(x){
  
  # List the distinct / unique values
  distinct_values <- unique(x)
  
  # Count the occurrence of each distinct value
  distinct_tabulate <- tabulate(match(x, distinct_values))
  
  # Return the value with the highest occurrence
  distinct_values[which.max(distinct_tabulate)]
}


train.set <- train.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))
test.set <- test.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))


#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

#mean imputation of mother age recruitment
train.set$mother_age_recruitment[is.na(train.set$mother_age_recruitment)] <- mean(train.set$mother_age_recruitment, na.rm = TRUE)
test.set$mother_age_recruitment[is.na(test.set$mother_age_recruitment)] <- mean(test.set$mother_age_recruitment, na.rm = TRUE)



### step 1
# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


### step 2
# for making sure train.set always has at least 2 samples for each category/ class for each column
# for columns that have a category/ class that has only 1 sample (eg X.83 only has 1 G G)
# add that category (eg G G) to that column in a new row
# and for the rest of the columns: add a random genotype already existing in the column

# Count the number of samples in each category class for each column
for (col in names(train.set)) {
  counts <- table(train.set[[col]])
  
  # Check if there is any category with only one sample
  if (sum(counts == 1) > 0) {
    # Print the column name if there is at least one category with one sample
    print(col)
  }
}

# Create a new row with all missing values
# ignore if no snps output from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 562))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 191] <- "A A" #190 + 1
new_row[, 192] <- "G G" #191 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(191,192)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:562, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(549), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



# Loop through the columns with 2 classes
# ignore if no output in the variable from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(test.set)))

for (col in cols_with_2_classes) {
  print(col)
  # Find the missing class in df1
  missing_class <- setdiff(levels(factor(train.set[[col]])), levels(factor(test.set[[col]])))
  
  
  #colnames(new_row) <- colnames(test.set)
  new_row[[col]] <- missing_class
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add the new row to test.set
# ignore if no output from prev code
colnames(new_row) <- colnames(test.set)
test.set <- rbind(test.set, new_row)

## for checking
# Make sure all columns in test.set have 3 unique values
for (col in names(test.set)) {
  if (length(unique(test.set[[col]])) < 3) {
    stop("Not all columns in test.set have 3 unique values.")
  }
}


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- seq (0, 1, by = 0.05) # 0.00 0.05 0.10 0.15 0.20 0.25 etc
lambda.vec <- 0:30 # 0  1  2  3  4  5  6 etc

# expand.grid: 
# Create a Data Frame from All Combinations of Factor Variables
elastic.grid <- expand.grid (alpha = alpha.vec, lambda = lambda.vec)
k = 10
repeats = 5

# repeatedcv:
# repeatedly performs X-fold cross-validation on the training data, 
# i.e. if you specify 5 repeats of 10-fold cross-validation, 
# it will perform 10-fold cross-validation on the training data 5 times, 
# using a different set of folds for each cross-validation
tr.Control <- trainControl (method = "repeatedcv",
                            number = k,
                            repeats = repeats,
                            search = "grid")


set.seed (1403)
elastic <- train (cg09322432 ~ ., data = train.set, 
                  method = 'glmnet', 
                  trControl = tr.Control, 
                  verbose = FALSE,
                  tuneGrid = elastic.grid)

elastic

elastic_trainresults <- data.frame(matrix(ncol = 4, nrow = 0))

elastic_trainresults <- cbind (elastic$results$alpha, elastic$results$lambda, elastic$results$RMSE, elastic$results$Rsquared)

colnames(elastic_trainresults) <- c('alpha', 'lambda', 'RMSE', 'Rsquared')

elastic_trainresults <- data.frame (elastic_trainresults)

# alpha = 0: ridge regression
# alpha = 1: lasso regression

plot (elastic)

lambda.opt <- elastic$bestTune$lambda
alpha.opt <- elastic$bestTune$alpha

index.opt <- which(elastic$results$lambda == lambda.opt  & 
                     elastic$results$alpha == alpha.opt )

results.elastic.opt <- elastic$results[index.opt, ]

results.elastic.opt


#predict on test.set without last row of dummy variables
elastic.pred <- predict(elastic, newdata = test.set[,1:561])


df <- data.frame(RMSE.el = RMSE(elastic.pred, test.set$cg09322432),
                 Rsquare.el = caret::R2(elastic.pred, test.set$cg09322432))



# XGBoost
#-----------------------------------------------------------------------------------------------
eval_results <- function (true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(RMSE = RMSE,
             Rsquare = R_square)
}
#-----------------------------------------------------------------------------------------------

# setting up xbg.DMatrix

# convert to sparse matrix
sparse_matrix <- sparse.model.matrix(cg09322432 ~ ., data = train.set)[,-1]
dtrain_mQTL = xgb.DMatrix(sparse_matrix, label = as.matrix(train.set$cg09322432))

### parameter tuning 100 rounds, random selection of parameters
best.param = list()
best.cv = Inf
best.cv.index = 0
best.nround = 0
set.seed(123)

for (iter in 1:50) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = sample(1:10, 1),
                eta = runif(1, 0.001, 0.05),
                subsample = runif(1, .5, .9),
                colsample_bytree = runif(1, .5, .9) 
  )
  cv.nround = 1000
  #uncomment nthread=16 if you want parallel processing, and you know the number of cores you have, i have 16 so i use 16 here
  mdcv <- xgb.cv(data=dtrain_mQTL, params = param, 
                 # nthread=6,
                 nfold=10, nrounds=cv.nround,
                 # verbose = FALSE, 
                 early_stopping_rounds=100)
  
  min.cv = min(mdcv$evaluation_log$test_rmse_mean)
  min.cv.index = which.min(mdcv$evaluation_log$test_rmse_mean)
  
  if (min.cv < best.cv) {
    best.cv = min.cv
    best.cv.index = min.cv.index
    best.param = param
  }
}

# best model 

best.param
nround = best.cv.index
nround

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)
xgb.save(best.boost, "xgboost.GTEx.model.non-genetic.probe3")


### predictions
xgboost_train_predictions <- as.data.frame(predict(best.boost, sparse_matrix))   

# convert test set to sparse matrix
sparse_matrix_test <- sparse.model.matrix(cg09322432 ~ ., data = test.set)[,-1]

xgboost_predictions <- as.data.frame(predict(best.boost, sparse_matrix_test))
test.set$cg09322432
xgboost_predictions

eval_results (train.set$cg09322432, xgboost_train_predictions, train.set) #train
eval_results (test.set$cg09322432, xgboost_predictions, test.set) #test

###################################
#                                 #
#      Probe 4 - cg10840135       #
#                                 #            
###################################
library(janitor) 
library(psych)
library(tidyverse)
library(xgboost)
library(DataExplorer)
library(caret)
library(randomForest)
library(readxl)
library(qdapRegex)
library(factoextra)
library(mgcv)
library(olsrr)
library(pls)
library(data.table)
library(leaps)
library(MASS)
library(neuralnet)
library(iml)
library(missForest)
library(Boruta)
library(glmnet)
library(ggplot2)
library(readxl)
library(dplyr)

library(readxl)

GenotypeFile <- read.csv("GTEx08snps_MAF10_genotype.csv")
pca <- read.csv("probe4_genotype_subset_LDP_pca.eigenvec_5PCs.csv")
demo <- read.csv("GUSTO_nongenetic_demo_withNA.csv")

Probe1 <- read.csv("Probe4_cg10840135.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")


#merge the files by FID
MainData <- GenotypeFile
MainData <- merge(MainData, demo, by = "FID")
MainData <- merge(MainData, pca, by = "FID")

MainData <- merge(MainData, Probe1, by = "FID")


test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#replace 0 with NA
train.set <- na_if(train.set, 0)
test.set <- na_if(test.set, 0)

# mode imputation
# highest no. of N/A in a col: 51
my_mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  if (length(xmode) > 1) xmode <- sample(xmode, 1)
  return(xmode)
}

for (var in 1:ncol(train.set)) {
  if (class(train.set[,var]) %in% c("character", "factor")) {
    train.set[is.na(train.set[,var]),var] <- my_mode(train.set[,var], na.rm = TRUE)
  }
}

for (var in 1:ncol(test.set)) {
  if (class(test.set[,var]) %in% c("character", "factor")) {
    test.set[is.na(test.set[,var]),var] <- my_mode(test.set[,var], na.rm = TRUE)
  }
}

# mode imputation for non-genetic cols
calc_mode <- function(x){
  
  # List the distinct / unique values
  distinct_values <- unique(x)
  
  # Count the occurrence of each distinct value
  distinct_tabulate <- tabulate(match(x, distinct_values))
  
  # Return the value with the highest occurrence
  distinct_values[which.max(distinct_tabulate)]
}


train.set <- train.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))
test.set <- test.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))


#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

#mean imputation of mother age recruitment
train.set$mother_age_recruitment[is.na(train.set$mother_age_recruitment)] <- mean(train.set$mother_age_recruitment, na.rm = TRUE)
test.set$mother_age_recruitment[is.na(test.set$mother_age_recruitment)] <- mean(test.set$mother_age_recruitment, na.rm = TRUE)



### step 1
# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


### step 2
# for making sure train.set always has at least 2 samples for each category/ class for each column
# for columns that have a category/ class that has only 1 sample (eg X.83 only has 1 G G)
# add that category (eg G G) to that column in a new row
# and for the rest of the columns: add a random genotype already existing in the column

# Count the number of samples in each category class for each column
for (col in names(train.set)) {
  counts <- table(train.set[[col]])
  
  # Check if there is any category with only one sample
  if (sum(counts == 1) > 0) {
    # Print the column name if there is at least one category with one sample
    print(col)
  }
}

# Create a new row with all missing values
# ignore if no snps output from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 562))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 191] <- "A A" #190 + 1
new_row[, 192] <- "G G" #191 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(191,192)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:562, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(549), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



# Loop through the columns with 2 classes
# ignore if no output in the variable from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(test.set)))

for (col in cols_with_2_classes) {
  print(col)
  # Find the missing class in df1
  missing_class <- setdiff(levels(factor(train.set[[col]])), levels(factor(test.set[[col]])))
  
  
  #colnames(new_row) <- colnames(test.set)
  new_row[[col]] <- missing_class
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add the new row to test.set
# ignore if no output from prev code
colnames(new_row) <- colnames(test.set)
test.set <- rbind(test.set, new_row)

## for checking
# Make sure all columns in test.set have 3 unique values
for (col in names(test.set)) {
  if (length(unique(test.set[[col]])) < 3) {
    stop("Not all columns in test.set have 3 unique values.")
  }
}


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- seq (0, 1, by = 0.05) # 0.00 0.05 0.10 0.15 0.20 0.25 etc
lambda.vec <- 0:30 # 0  1  2  3  4  5  6 etc

# expand.grid: 
# Create a Data Frame from All Combinations of Factor Variables
elastic.grid <- expand.grid (alpha = alpha.vec, lambda = lambda.vec)
k = 10
repeats = 5

# repeatedcv:
# repeatedly performs X-fold cross-validation on the training data, 
# i.e. if you specify 5 repeats of 10-fold cross-validation, 
# it will perform 10-fold cross-validation on the training data 5 times, 
# using a different set of folds for each cross-validation
tr.Control <- trainControl (method = "repeatedcv",
                            number = k,
                            repeats = repeats,
                            search = "grid")


set.seed (1403)
elastic <- train (cg10840135 ~ ., data = train.set, 
                  method = 'glmnet', 
                  trControl = tr.Control, 
                  verbose = FALSE,
                  tuneGrid = elastic.grid)

elastic

elastic_trainresults <- data.frame(matrix(ncol = 4, nrow = 0))

elastic_trainresults <- cbind (elastic$results$alpha, elastic$results$lambda, elastic$results$RMSE, elastic$results$Rsquared)

colnames(elastic_trainresults) <- c('alpha', 'lambda', 'RMSE', 'Rsquared')

elastic_trainresults <- data.frame (elastic_trainresults)

# alpha = 0: ridge regression
# alpha = 1: lasso regression

plot (elastic)

lambda.opt <- elastic$bestTune$lambda
alpha.opt <- elastic$bestTune$alpha

index.opt <- which(elastic$results$lambda == lambda.opt  & 
                     elastic$results$alpha == alpha.opt )

results.elastic.opt <- elastic$results[index.opt, ]

results.elastic.opt


#predict on test.set without last row of dummy variables
elastic.pred <- predict(elastic, newdata = test.set[,1:561])


df <- data.frame(RMSE.el = RMSE(elastic.pred, test.set$cg10840135),
                 Rsquare.el = caret::R2(elastic.pred, test.set$cg10840135))



# XGBoost
#-----------------------------------------------------------------------------------------------
eval_results <- function (true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(RMSE = RMSE,
             Rsquare = R_square)
}
#-----------------------------------------------------------------------------------------------

# setting up xbg.DMatrix

# convert to sparse matrix
sparse_matrix <- sparse.model.matrix(cg10840135 ~ ., data = train.set)[,-1]
dtrain_mQTL = xgb.DMatrix(sparse_matrix, label = as.matrix(train.set$cg10840135))

### parameter tuning 100 rounds, random selection of parameters
best.param = list()
best.cv = Inf
best.cv.index = 0
best.nround = 0
set.seed(123)

for (iter in 1:50) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = sample(1:10, 1),
                eta = runif(1, 0.001, 0.05),
                subsample = runif(1, .5, .9),
                colsample_bytree = runif(1, .5, .9) 
  )
  cv.nround = 1000
  #uncomment nthread=16 if you want parallel processing, and you know the number of cores you have, i have 16 so i use 16 here
  mdcv <- xgb.cv(data=dtrain_mQTL, params = param, 
                 # nthread=6,
                 nfold=10, nrounds=cv.nround,
                 # verbose = FALSE, 
                 early_stopping_rounds=100)
  
  min.cv = min(mdcv$evaluation_log$test_rmse_mean)
  min.cv.index = which.min(mdcv$evaluation_log$test_rmse_mean)
  
  if (min.cv < best.cv) {
    best.cv = min.cv
    best.cv.index = min.cv.index
    best.param = param
  }
}

# best model 

best.param
nround = best.cv.index
nround

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)
xgb.save(best.boost, "xgboost.GTEx.model.non-genetic.probe4")


### predictions
xgboost_train_predictions <- as.data.frame(predict(best.boost, sparse_matrix))   

# convert test set to sparse matrix
sparse_matrix_test <- sparse.model.matrix(cg10840135 ~ ., data = test.set)[,-1]

xgboost_predictions <- as.data.frame(predict(best.boost, sparse_matrix_test))
test.set$cg10840135
xgboost_predictions

eval_results (train.set$cg10840135, xgboost_train_predictions, train.set) #train
eval_results (test.set$cg10840135, xgboost_predictions, test.set) #test

###################################
#                                 #
#      Probe 5 - cg15597984       #
#                                 #            
###################################
library(janitor) 
library(psych)
library(tidyverse)
library(xgboost)
library(DataExplorer)
library(caret)
library(randomForest)
library(readxl)
library(qdapRegex)
library(factoextra)
library(mgcv)
library(olsrr)
library(pls)
library(data.table)
library(leaps)
library(MASS)
library(neuralnet)
library(iml)
library(missForest)
library(Boruta)
library(glmnet)
library(ggplot2)
library(readxl)
library(dplyr)

library(readxl)

GenotypeFile <- read.csv("GTEx08snps_MAF10_genotype.csv")
pca <- read.csv("probe5_genotype_subset_LDP_pca.eigenvec_5PCs.csv")
demo <- read.csv("GUSTO_nongenetic_demo_withNA.csv")

Probe1 <- read.csv("Probe5_cg15597984.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")


#merge the files by FID
MainData <- GenotypeFile
MainData <- merge(MainData, demo, by = "FID")
MainData <- merge(MainData, pca, by = "FID")

MainData <- merge(MainData, Probe1, by = "FID")


test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#replace 0 with NA
train.set <- na_if(train.set, 0)
test.set <- na_if(test.set, 0)

# mode imputation
# highest no. of N/A in a col: 51
my_mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  if (length(xmode) > 1) xmode <- sample(xmode, 1)
  return(xmode)
}

for (var in 1:ncol(train.set)) {
  if (class(train.set[,var]) %in% c("character", "factor")) {
    train.set[is.na(train.set[,var]),var] <- my_mode(train.set[,var], na.rm = TRUE)
  }
}

for (var in 1:ncol(test.set)) {
  if (class(test.set[,var]) %in% c("character", "factor")) {
    test.set[is.na(test.set[,var]),var] <- my_mode(test.set[,var], na.rm = TRUE)
  }
}

# mode imputation for non-genetic cols
calc_mode <- function(x){
  
  # List the distinct / unique values
  distinct_values <- unique(x)
  
  # Count the occurrence of each distinct value
  distinct_tabulate <- tabulate(match(x, distinct_values))
  
  # Return the value with the highest occurrence
  distinct_values[which.max(distinct_tabulate)]
}


train.set <- train.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))
test.set <- test.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))


#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

#mean imputation of mother age recruitment
train.set$mother_age_recruitment[is.na(train.set$mother_age_recruitment)] <- mean(train.set$mother_age_recruitment, na.rm = TRUE)
test.set$mother_age_recruitment[is.na(test.set$mother_age_recruitment)] <- mean(test.set$mother_age_recruitment, na.rm = TRUE)



### step 1
# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


### step 2
# for making sure train.set always has at least 2 samples for each category/ class for each column
# for columns that have a category/ class that has only 1 sample (eg X.83 only has 1 G G)
# add that category (eg G G) to that column in a new row
# and for the rest of the columns: add a random genotype already existing in the column

# Count the number of samples in each category class for each column
for (col in names(train.set)) {
  counts <- table(train.set[[col]])
  
  # Check if there is any category with only one sample
  if (sum(counts == 1) > 0) {
    # Print the column name if there is at least one category with one sample
    print(col)
  }
}

# Create a new row with all missing values
# ignore if no snps output from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 562))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 191] <- "A A" #190 + 1
new_row[, 192] <- "G G" #191 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(191,192)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:562, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(549), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



# Loop through the columns with 2 classes
# ignore if no output in the variable from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(test.set)))

for (col in cols_with_2_classes) {
  print(col)
  # Find the missing class in df1
  missing_class <- setdiff(levels(factor(train.set[[col]])), levels(factor(test.set[[col]])))
  
  
  #colnames(new_row) <- colnames(test.set)
  new_row[[col]] <- missing_class
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add the new row to test.set
# ignore if no output from prev code
colnames(new_row) <- colnames(test.set)
test.set <- rbind(test.set, new_row)

## for checking
# Make sure all columns in test.set have 3 unique values
for (col in names(test.set)) {
  if (length(unique(test.set[[col]])) < 3) {
    stop("Not all columns in test.set have 3 unique values.")
  }
}


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- seq (0, 1, by = 0.05) # 0.00 0.05 0.10 0.15 0.20 0.25 etc
lambda.vec <- 0:30 # 0  1  2  3  4  5  6 etc

# expand.grid: 
# Create a Data Frame from All Combinations of Factor Variables
elastic.grid <- expand.grid (alpha = alpha.vec, lambda = lambda.vec)
k = 10
repeats = 5

# repeatedcv:
# repeatedly performs X-fold cross-validation on the training data, 
# i.e. if you specify 5 repeats of 10-fold cross-validation, 
# it will perform 10-fold cross-validation on the training data 5 times, 
# using a different set of folds for each cross-validation
tr.Control <- trainControl (method = "repeatedcv",
                            number = k,
                            repeats = repeats,
                            search = "grid")


set.seed (1403)
elastic <- train (cg15597984 ~ ., data = train.set, 
                  method = 'glmnet', 
                  trControl = tr.Control, 
                  verbose = FALSE,
                  tuneGrid = elastic.grid)

elastic

elastic_trainresults <- data.frame(matrix(ncol = 4, nrow = 0))

elastic_trainresults <- cbind (elastic$results$alpha, elastic$results$lambda, elastic$results$RMSE, elastic$results$Rsquared)

colnames(elastic_trainresults) <- c('alpha', 'lambda', 'RMSE', 'Rsquared')

elastic_trainresults <- data.frame (elastic_trainresults)

# alpha = 0: ridge regression
# alpha = 1: lasso regression

plot (elastic)

lambda.opt <- elastic$bestTune$lambda
alpha.opt <- elastic$bestTune$alpha

index.opt <- which(elastic$results$lambda == lambda.opt  & 
                     elastic$results$alpha == alpha.opt )

results.elastic.opt <- elastic$results[index.opt, ]

results.elastic.opt


#predict on test.set without last row of dummy variables
elastic.pred <- predict(elastic, newdata = test.set[,1:561])


df <- data.frame(RMSE.el = RMSE(elastic.pred, test.set$cg15597984),
                 Rsquare.el = caret::R2(elastic.pred, test.set$cg15597984))



# XGBoost
#-----------------------------------------------------------------------------------------------
eval_results <- function (true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(RMSE = RMSE,
             Rsquare = R_square)
}
#-----------------------------------------------------------------------------------------------

# setting up xbg.DMatrix

# convert to sparse matrix
sparse_matrix <- sparse.model.matrix(cg15597984 ~ ., data = train.set)[,-1]
dtrain_mQTL = xgb.DMatrix(sparse_matrix, label = as.matrix(train.set$cg15597984))

### parameter tuning 100 rounds, random selection of parameters
best.param = list()
best.cv = Inf
best.cv.index = 0
best.nround = 0
set.seed(123)

for (iter in 1:50) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = sample(1:10, 1),
                eta = runif(1, 0.001, 0.05),
                subsample = runif(1, .5, .9),
                colsample_bytree = runif(1, .5, .9) 
  )
  cv.nround = 1000
  #uncomment nthread=16 if you want parallel processing, and you know the number of cores you have, i have 16 so i use 16 here
  mdcv <- xgb.cv(data=dtrain_mQTL, params = param, 
                 # nthread=6,
                 nfold=10, nrounds=cv.nround,
                 # verbose = FALSE, 
                 early_stopping_rounds=100)
  
  min.cv = min(mdcv$evaluation_log$test_rmse_mean)
  min.cv.index = which.min(mdcv$evaluation_log$test_rmse_mean)
  
  if (min.cv < best.cv) {
    best.cv = min.cv
    best.cv.index = min.cv.index
    best.param = param
  }
}

# best model 

best.param
nround = best.cv.index
nround

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)
xgb.save(best.boost, "xgboost.GTEx.model.non-genetic.probe5")


### predictions
xgboost_train_predictions <- as.data.frame(predict(best.boost, sparse_matrix))   

# convert test set to sparse matrix
sparse_matrix_test <- sparse.model.matrix(cg15597984 ~ ., data = test.set)[,-1]

xgboost_predictions <- as.data.frame(predict(best.boost, sparse_matrix_test))
test.set$cg15597984
xgboost_predictions

eval_results (train.set$cg15597984, xgboost_train_predictions, train.set) #train
eval_results (test.set$cg15597984, xgboost_predictions, test.set) #test

###################################
#                                 #
#      Probe 7 - cg20046859       #
#                                 #            
###################################
library(janitor) 
library(psych)
library(tidyverse)
library(xgboost)
library(DataExplorer)
library(caret)
library(randomForest)
library(readxl)
library(qdapRegex)
library(factoextra)
library(mgcv)
library(olsrr)
library(pls)
library(data.table)
library(leaps)
library(MASS)
library(neuralnet)
library(iml)
library(missForest)
library(Boruta)
library(glmnet)
library(ggplot2)
library(readxl)
library(dplyr)

library(readxl)

GenotypeFile <- read.csv("GTEx08snps_MAF10_genotype.csv")
pca <- read.csv("probe7_genotype_subset_LDP_pca.eigenvec_5PCs.csv")
demo <- read.csv("GUSTO_nongenetic_demo_withNA.csv")

Probe1 <- read.csv("Probe7_cg20046859.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")


#merge the files by FID
MainData <- GenotypeFile
MainData <- merge(MainData, demo, by = "FID")
MainData <- merge(MainData, pca, by = "FID")

MainData <- merge(MainData, Probe1, by = "FID")


test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#replace 0 with NA
train.set <- na_if(train.set, 0)
test.set <- na_if(test.set, 0)

# mode imputation
# highest no. of N/A in a col: 51
my_mode <- function (x, na.rm) {
  xtab <- table(x)
  xmode <- names(which(xtab == max(xtab)))
  if (length(xmode) > 1) xmode <- sample(xmode, 1)
  return(xmode)
}

for (var in 1:ncol(train.set)) {
  if (class(train.set[,var]) %in% c("character", "factor")) {
    train.set[is.na(train.set[,var]),var] <- my_mode(train.set[,var], na.rm = TRUE)
  }
}

for (var in 1:ncol(test.set)) {
  if (class(test.set[,var]) %in% c("character", "factor")) {
    test.set[is.na(test.set[,var]),var] <- my_mode(test.set[,var], na.rm = TRUE)
  }
}

# mode imputation for non-genetic cols
calc_mode <- function(x){
  
  # List the distinct / unique values
  distinct_values <- unique(x)
  
  # Count the occurrence of each distinct value
  distinct_tabulate <- tabulate(match(x, distinct_values))
  
  # Return the value with the highest occurrence
  distinct_values[which.max(distinct_tabulate)]
}


train.set <- train.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))
test.set <- test.set %>%  mutate(across(everything(), ~replace_na(.x, calc_mode(.x))))


#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

#mean imputation of mother age recruitment
train.set$mother_age_recruitment[is.na(train.set$mother_age_recruitment)] <- mean(train.set$mother_age_recruitment, na.rm = TRUE)
test.set$mother_age_recruitment[is.na(test.set$mother_age_recruitment)] <- mean(test.set$mother_age_recruitment, na.rm = TRUE)



### step 1
# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


### step 2
# for making sure train.set always has at least 2 samples for each category/ class for each column
# for columns that have a category/ class that has only 1 sample (eg X.83 only has 1 G G)
# add that category (eg G G) to that column in a new row
# and for the rest of the columns: add a random genotype already existing in the column

# Count the number of samples in each category class for each column
for (col in names(train.set)) {
  counts <- table(train.set[[col]])
  
  # Check if there is any category with only one sample
  if (sum(counts == 1) > 0) {
    # Print the column name if there is at least one category with one sample
    print(col)
  }
}

# Create a new row with all missing values
# ignore if no snps output from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 562))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 191] <- "A A" #190 + 1
new_row[, 192] <- "G G" #191 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(191,192)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:562, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(549), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



# Loop through the columns with 2 classes
# ignore if no output in the variable from prev code
new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(test.set)))

for (col in cols_with_2_classes) {
  print(col)
  # Find the missing class in df1
  missing_class <- setdiff(levels(factor(train.set[[col]])), levels(factor(test.set[[col]])))
  
  
  #colnames(new_row) <- colnames(test.set)
  new_row[[col]] <- missing_class
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add random values to the rest of the columns in the new row
# ignore if no output from prev code
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add the new row to test.set
# ignore if no output from prev code
colnames(new_row) <- colnames(test.set)
test.set <- rbind(test.set, new_row)

## for checking
# Make sure all columns in test.set have 3 unique values
for (col in names(test.set)) {
  if (length(unique(test.set[[col]])) < 3) {
    stop("Not all columns in test.set have 3 unique values.")
  }
}


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- seq (0, 1, by = 0.05) # 0.00 0.05 0.10 0.15 0.20 0.25 etc
lambda.vec <- 0:30 # 0  1  2  3  4  5  6 etc

# expand.grid: 
# Create a Data Frame from All Combinations of Factor Variables
elastic.grid <- expand.grid (alpha = alpha.vec, lambda = lambda.vec)
k = 10
repeats = 5

# repeatedcv:
# repeatedly performs X-fold cross-validation on the training data, 
# i.e. if you specify 5 repeats of 10-fold cross-validation, 
# it will perform 10-fold cross-validation on the training data 5 times, 
# using a different set of folds for each cross-validation
tr.Control <- trainControl (method = "repeatedcv",
                            number = k,
                            repeats = repeats,
                            search = "grid")


set.seed (1403)
elastic <- train (cg20046859 ~ ., data = train.set, 
                  method = 'glmnet', 
                  trControl = tr.Control, 
                  verbose = FALSE,
                  tuneGrid = elastic.grid)

elastic

elastic_trainresults <- data.frame(matrix(ncol = 4, nrow = 0))

elastic_trainresults <- cbind (elastic$results$alpha, elastic$results$lambda, elastic$results$RMSE, elastic$results$Rsquared)

colnames(elastic_trainresults) <- c('alpha', 'lambda', 'RMSE', 'Rsquared')

elastic_trainresults <- data.frame (elastic_trainresults)

# alpha = 0: ridge regression
# alpha = 1: lasso regression

plot (elastic)

lambda.opt <- elastic$bestTune$lambda
alpha.opt <- elastic$bestTune$alpha

index.opt <- which(elastic$results$lambda == lambda.opt  & 
                     elastic$results$alpha == alpha.opt )

results.elastic.opt <- elastic$results[index.opt, ]

results.elastic.opt


#predict on test.set without last row of dummy variables
elastic.pred <- predict(elastic, newdata = test.set[,1:561])


df <- data.frame(RMSE.el = RMSE(elastic.pred, test.set$cg20046859),
                 Rsquare.el = caret::R2(elastic.pred, test.set$cg20046859))



# XGBoost
#-----------------------------------------------------------------------------------------------
eval_results <- function (true, predicted, df) {
  SSE <- sum((predicted - true)^2)
  SST <- sum((true - mean(true))^2)
  R_square <- 1 - SSE / SST
  RMSE = sqrt(SSE/nrow(df))
  # Model performance metrics
  data.frame(RMSE = RMSE,
             Rsquare = R_square)
}
#-----------------------------------------------------------------------------------------------

# setting up xbg.DMatrix

# convert to sparse matrix
sparse_matrix <- sparse.model.matrix(cg20046859 ~ ., data = train.set)[,-1]
dtrain_mQTL = xgb.DMatrix(sparse_matrix, label = as.matrix(train.set$cg20046859))

### parameter tuning 100 rounds, random selection of parameters
best.param = list()
best.cv = Inf
best.cv.index = 0
best.nround = 0
set.seed(123)

for (iter in 1:50) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = sample(1:10, 1),
                eta = runif(1, 0.001, 0.05),
                subsample = runif(1, .5, .9),
                colsample_bytree = runif(1, .5, .9) 
  )
  cv.nround = 1000
  #uncomment nthread=16 if you want parallel processing, and you know the number of cores you have, i have 16 so i use 16 here
  mdcv <- xgb.cv(data=dtrain_mQTL, params = param, 
                 # nthread=6,
                 nfold=10, nrounds=cv.nround,
                 # verbose = FALSE, 
                 early_stopping_rounds=100)
  
  min.cv = min(mdcv$evaluation_log$test_rmse_mean)
  min.cv.index = which.min(mdcv$evaluation_log$test_rmse_mean)
  
  if (min.cv < best.cv) {
    best.cv = min.cv
    best.cv.index = min.cv.index
    best.param = param
  }
}

# best model 

best.param
nround = best.cv.index
nround

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)
xgb.save(best.boost, "xgboost.GTEx.model.non-genetic.probe7")


### predictions
xgboost_train_predictions <- as.data.frame(predict(best.boost, sparse_matrix))   

# convert test set to sparse matrix
sparse_matrix_test <- sparse.model.matrix(cg20046859 ~ ., data = test.set)[,-1]

xgboost_predictions <- as.data.frame(predict(best.boost, sparse_matrix_test))
test.set$cg20046859
xgboost_predictions

eval_results (train.set$cg20046859, xgboost_train_predictions, train.set) #train
eval_results (test.set$cg20046859, xgboost_predictions, test.set) #test

