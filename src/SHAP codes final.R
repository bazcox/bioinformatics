## SHAP 

###################################
#                                 #
#      Probe 1 - cg04692870       #
#                                 #            
###################################

# best model: Elastic net Ind GWAS after B-H genetics only (2 SNPs) &
#             XGBoost Ind GWAS after B-H genetics only (2 SNPs)

GenotypeFile <- read.csv("GWAS_trainset_afterBH_indiv_probe1_v2_genotype.csv",na = c("", "NA", "N/A",0))

colnames(GenotypeFile)[ncol(GenotypeFile) - 1] <- "rs133335"
colnames(GenotypeFile)[ncol(GenotypeFile)] <- "rs133344"


Probe1 <- read.csv("Probe1_cg04692870_v2.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")

#merge the files by FID
MainData <- merge(GenotypeFile, Probe1, by = "FID")

test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

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

# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


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


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(test.set, function(x) length(unique(x))) == 2)

write.csv(train.set, "Probe1_top_snps_EL_XGBoost_trainset.csv")

#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- 0
lambda.vec <- 0

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

################# DALEX ########################

library(DALEX)

explainer <- explain(elastic, data = train.set[, -which(names(train.set) == "cg04692870")], y = train.set$cg04692870)

# Calculate feature importance using variable importance profiles (VIP)
vip_data <- variable_importance(explainer)

plot_data <- vip_data[, c("variable", "dropout_loss")]
plot_data <- plot_data[order(-plot_data$dropout_loss), ]

ggplot(plot_data, aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Dropout Loss", title = "Variable Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()


############################# SHAP for elastic net ########################
library(ggplot2)

set.seed(309)
elastic_coef <- coef(elastic$finalModel, s = 0) # change s = 0 to your best lambda value

# Assuming 'elastic_coef' is the coefficient matrix obtained from the Elastic Net model
feature_importance <- data.frame(
  Feature = rownames(elastic_coef),
  Coefficient = as.numeric(elastic_coef)
)

# Take absolute values of the coefficients for importance
feature_importance$AbsoluteCoefficient <- abs(feature_importance$Coefficient)

# Remove the intercept (if it exists)
feature_importance <- feature_importance[feature_importance$Feature != "(Intercept)",]

# Sort the features by their absolute coefficients in descending order
feature_importance <- feature_importance[order(-feature_importance$AbsoluteCoefficient),]

# Select the top 10 features
top_10_features <- head(feature_importance, 10)

# Create a bar plot for the top 10 features
ggplot(top_10_features, aes(x = reorder(Feature, AbsoluteCoefficient), y = AbsoluteCoefficient)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Feature", y = "Importance") +
  theme_minimal()


# XGBoost

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

for (iter in 1:1) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = 10,
                eta = 0.003009989,
                subsample = 0.8383967,
                colsample_bytree = 0.5493245
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
nround = 1000

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)

################ SHAP ######################

set.seed(309)
shap.score.rank <- function(xgb_model = xgb_mod, shap_approx = TRUE, 
                            X_train = mydata$train_mm){
  require(xgboost)
  require(data.table)
  shap_contrib <- predict(xgb_model, X_train,
                          predcontrib = TRUE, approxcontrib = shap_approx)
  shap_contrib <- as.data.table(shap_contrib)
  shap_contrib[,BIAS:=NULL]
  cat('make SHAP score by decreasing order\n\n')
  mean_shap_score <- colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = T)]
  return(list(shap_score = shap_contrib,
              mean_shap_score = (mean_shap_score)))
}

plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce') # for `geom_sina`
  plot1 <- ggplot(data = data_long)+
    coord_flip() + 
    # sina plot: 
    geom_sina(aes(x = variable, y = value, color = stdfvalue)) +
    # print the mean absolute value: 
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) + 
    scale_color_gradient(low="#FFCC33", high="#6600CC", 
                         breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), # remove axis line
          legend.position="bottom") + 
    geom_hline(yintercept = 0) + # the vertical line
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}

# a function to standardize feature values into same range
std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}


# prep shap data
shap.prep <- function(shap  = shap_result, X_train = mydata$train_mm, top_n){
  require(ggforce)
  # descending order
  if (missing(top_n)) top_n <- dim(X_train)[2] # by default, use all features
  if (!top_n%in%c(1:dim(X_train)[2])) stop('supply correct top_n')
  require(data.table)
  shap_score_sub <- as.data.table(shap$shap_score)
  shap_score_sub <- shap_score_sub[, names(shap$mean_shap_score)[1:top_n], with = F]
  shap_score_long <- melt.data.table(shap_score_sub, measure.vars = colnames(shap_score_sub))
  
  # feature values: the values in the original dataset
  fv_sub <- as.data.table(X_train)[, names(shap$mean_shap_score)[1:top_n], with = F]
  # standardize feature values
  fv_sub_long <- melt.data.table(fv_sub, measure.vars = colnames(fv_sub))
  fv_sub_long[, stdfvalue := std1(value), by = "variable"]
  # SHAP value: value
  # raw feature value: rfvalue; 
  # standarized: stdfvalue
  names(fv_sub_long) <- c("variable", "rfvalue", "stdfvalue" )
  shap_long2 <- cbind(shap_score_long, fv_sub_long[,c('rfvalue','stdfvalue')])
  shap_long2[, mean_value := mean(abs(value)), by = variable]
  setkey(shap_long2, variable)
  return(shap_long2) 
}

var_importance <- function(shap_result, top_n = 10) {
  # Extract feature names and SHAP values
  feature_names <- names(shap_result$mean_shap_score)
  shap_values <- shap_result$mean_shap_score
  
  # Create a data frame with feature names and SHAP values
  shap_df <- tibble(var = feature_names, importance = shap_values)
  
  # Remove rows with NA values
  shap_df <- shap_df[complete.cases(shap_df), ]
  
  # Take the top_n features
  top_features <- shap_df %>%
    arrange(desc(abs(importance))) %>%
    head(top_n)
  
  return(top_features)
}

# Modified function to return top 10 SHAP values as a data frame
top_n_shap_df <- function(shap_result, top_n = 10) {
  feature_names <- names(shap_result$mean_shap_score)
  shap_values <- shap_result$mean_shap_score
  shap_df <- tibble(var = feature_names, importance = shap_values)
  
  # Remove rows with NA values
  shap_df <- shap_df[complete.cases(shap_df), ]
  
  # Take the top_n features
  top_features <- shap_df %>%
    arrange(desc(abs(importance))) %>%
    head(top_n)
  
  return(top_features)
}

#shap_values=predict(best.boost, dtrain_mQTL, predcontrib = TRUE, approxcontrib = F)

shap_result_bike = shap.score.rank(xgb_model = best.boost, 
                                   X_train =dtrain_mQTL,
                                   shap_approx = F
)

# Get the top 10 SHAP values as a data frame
top_10_shap_df <- top_n_shap_df(shap_result_bike, top_n = 10)
print(top_10_shap_df)

# Plot var importance based on SHAP
ggplot(top_10_shap_df, aes(x = reorder(var, importance), y = importance)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_light() + 
  theme(axis.title.y = element_blank()) +
  labs(y = "Importance", x = "Feature")

##################################



###################################
#                                 #
#      Probe 3 - cg09322432       #
#                                 #            
###################################

# best model: Elastic net Ind GWAS uncorrected genetics + non-genetics

GenotypeFile <- read.csv("GWAS_trainset_beforeBH_indiv_probe3_snps_MAF10_genotype.csv",na = c("", "NA", "N/A",0))
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
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 204))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 53] <- "C C" #52 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(53)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:204, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(191), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



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

write.csv(train.set,"Probe3_top_snps_trainset.csv", row.names = FALSE)

#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- 0
lambda.vec <- 0

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


################# DALEX ########################

library(DALEX)

explainer <- explain(elastic, data = train.set[, -which(names(train.set) == "cg09322432")], y = train.set$cg09322432)

# Calculate feature importance using variable importance profiles (VIP)
vip_data <- variable_importance(explainer)

plot_data <- vip_data[, c("variable", "dropout_loss")]
plot_data <- plot_data[order(-plot_data$dropout_loss), ]

top_40_plot_data <- head(plot_data, 40)

ggplot(top_40_plot_data, aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Dropout Loss", title = "Top 10 Variable Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

############################# SHAP for elastic net ########################
library(ggplot2)

set.seed(1403)
elastic_coef <- coef(elastic$finalModel, s = 0) # change s = 0 to your best lambda value

# Assuming 'elastic_coef' is the coefficient matrix obtained from the Elastic Net model
feature_importance <- data.frame(
  Feature = rownames(elastic_coef),
  Coefficient = as.numeric(elastic_coef)
)

# Take absolute values of the coefficients for importance
feature_importance$AbsoluteCoefficient <- abs(feature_importance$Coefficient)

# Remove the intercept (if it exists)
feature_importance <- feature_importance[feature_importance$Feature != "(Intercept)",]

# Sort the features by their absolute coefficients in descending order
feature_importance <- feature_importance[order(-feature_importance$AbsoluteCoefficient),]

# Select the top 10 features
top_10_features <- head(feature_importance, 10)

# Create a bar plot for the top 10 features
ggplot(top_10_features, aes(x = reorder(Feature, AbsoluteCoefficient), y = AbsoluteCoefficient)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Feature", y = "Importance") +
  theme_minimal()



###################################
#                                 #
#      Probe 4 - cg10840135       #
#                                 #            
###################################

# best model: Elastic Nets GTEx (548 snps) &
#             XGBoost 2MB genetic only (2406 snps)

GenotypeFile <- read.csv("GTEx08snps_MAF10_genotype.csv",na = c("", "NA", "N/A",0))
Probe1 <- read.csv("Probe4_cg10840135.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")

#merge the files by FID
MainData <- merge(GenotypeFile, Probe1, by = "FID")

test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]


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



## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(test.set, function(x) length(unique(x))) == 2)

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

write.csv(train.set,"Probe4_top_snps_trainset.csv", row.names = FALSE)

#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- 0
lambda.vec <- 1

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

################# DALEX ########################

library(DALEX)

explainer <- explain(elastic, data = train.set[, -which(names(train.set) == "cg10840135")], y = train.set$cg10840135)

# Calculate feature importance using variable importance profiles (VIP)
vip_data <- variable_importance(explainer)

plot_data <- vip_data[, c("variable", "dropout_loss")]
plot_data <- plot_data[order(-plot_data$dropout_loss), ]

top_40_plot_data <- head(plot_data, 39)

ggplot(top_40_plot_data, aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Dropout Loss", title = "Top 10 Variable Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

############################# SHAP for elastic net ########################
library(ggplot2)

elastic_coef <- coef(elastic$finalModel, s = 1) # change s = 0 to your best lambda value

# Assuming 'elastic_coef' is the coefficient matrix obtained from the Elastic Net model
feature_importance <- data.frame(
  Feature = rownames(elastic_coef),
  Coefficient = as.numeric(elastic_coef)
)

# Take absolute values of the coefficients for importance
feature_importance$AbsoluteCoefficient <- abs(feature_importance$Coefficient)

# Remove the intercept (if it exists)
feature_importance <- feature_importance[feature_importance$Feature != "(Intercept)",]

# Sort the features by their absolute coefficients in descending order
feature_importance <- feature_importance[order(-feature_importance$AbsoluteCoefficient),]

# Select the top 10 features
top_10_features <- head(feature_importance, 10)

# Create a bar plot for the top 10 features
ggplot(top_10_features, aes(x = reorder(Feature, AbsoluteCoefficient), y = AbsoluteCoefficient)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Feature", y = "Importance") +
  theme_minimal()


## XGBoost
GenotypeFile <- read.csv("08-41.5MB-43.6MBsnps_MAF10_genotype.csv",na = c("", "NA", "N/A",0))
Probe1 <- read.csv("Probe4_cg10840135.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")

#merge the files by FID
MainData <- merge(GenotypeFile, Probe1, by = "FID")

test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]


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


# Define the columns to exclude
# ignore if no output from prev code & go to step 2
exclude_cols <- c("rs9607897","rs9607991")

# Initialize a flag to indicate whether a row has been moved
# ignore if no output from prev code
row_moved <- FALSE

# Check each excluded column for "AA"
# ignore if no output from prev code
for (col in exclude_cols) {
  # Initialize a variable to store the desired value for the current column
  desired_value <- ""
  
  # Determine the desired value based on the current column
  if (col == "rs9607897") {
    desired_value <- "T T"
  } else if (col == "rs9607991") {
    desired_value <- "T T"
  }
  
  # Find rows in test.set that match the desired value
  rows_to_move <- which(test.set[, col] == desired_value)
  
  # If there are rows with the desired value in the current column
  if (length(rows_to_move) > 0) {
    # If there's only one row to move, use it
    if (length(rows_to_move) == 1) {
      row_index <- rows_to_move
    } else {
      # Select a random row
      row_index <- sample(rows_to_move, 1)
    }
    
    # Add the selected row to train.set
    train.set <- rbind(train.set, test.set[row_index, ])
    
    # Remove the selected row from test.set
    test.set <- test.set[-row_index, ]
    
    # Set the flag to indicate that a row has been moved
    row_moved <- TRUE
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
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 2407))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 
# change the letters too! 
# ignore if no output from prev code

#                     actual column + 1
new_row[, 1178] <- "T T" #1177 + 1
new_row[, 1242] <- "C C" #1241 + 1
new_row[, 2148] <- "A A" #2147 + 1
new_row[, 2170] <- "A A" #2169 + 1
new_row[, 2183] <- "T T" #2182 + 1
new_row[, 2188] <- "T T" #2187 + 1
new_row[, 2406] <- "C C" #2405 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(1178,1242,2148,2170,2183,2188,2406)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:2407, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
# ignore if no output from prev code
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(test.set, function(x) length(unique(x))) == 2)

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


write.csv(train.set,"Probe4_XGBoost_top_snps_trainset.csv", row.names = FALSE)


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

for (iter in 1:1) {
  param <- list(objective = "reg:squarederror", ###### this is squared error. 
                max_depth = 1,
                eta = 0.006558051,
                subsample = 0.7592239,
                colsample_bytree = 0.7361971 
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
nround = 924
nround

set.seed(123)
best.boost <- xgboost(data=dtrain_mQTL, params=best.param, nrounds=nround, verbose = FALSE)

################ SHAP ######################
shap.score.rank <- function(xgb_model = xgb_mod, shap_approx = TRUE, 
                            X_train = mydata$train_mm){
  require(xgboost)
  require(data.table)
  shap_contrib <- predict(xgb_model, X_train,
                          predcontrib = TRUE, approxcontrib = shap_approx)
  shap_contrib <- as.data.table(shap_contrib)
  shap_contrib[,BIAS:=NULL]
  cat('make SHAP score by decreasing order\n\n')
  mean_shap_score <- colMeans(abs(shap_contrib))[order(colMeans(abs(shap_contrib)), decreasing = T)]
  return(list(shap_score = shap_contrib,
              mean_shap_score = (mean_shap_score)))
}

plot.shap.summary <- function(data_long){
  x_bound <- max(abs(data_long$value))
  require('ggforce') # for `geom_sina`
  plot1 <- ggplot(data = data_long)+
    coord_flip() + 
    # sina plot: 
    geom_sina(aes(x = variable, y = value, color = stdfvalue)) +
    # print the mean absolute value: 
    geom_text(data = unique(data_long[, c("variable", "mean_value"), with = F]),
              aes(x = variable, y=-Inf, label = sprintf("%.3f", mean_value)),
              size = 3, alpha = 0.7,
              hjust = -0.2, 
              fontface = "bold") + # bold
    # # add a "SHAP" bar notation
    # annotate("text", x = -Inf, y = -Inf, vjust = -0.2, hjust = 0, size = 3,
    #          label = expression(group("|", bar(SHAP), "|"))) + 
    scale_color_gradient(low="#FFCC33", high="#6600CC", 
                         breaks=c(0,1), labels=c("Low","High")) +
    theme_bw() + 
    theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), # remove axis line
          legend.position="bottom") + 
    geom_hline(yintercept = 0) + # the vertical line
    scale_y_continuous(limits = c(-x_bound, x_bound)) +
    # reverse the order of features
    scale_x_discrete(limits = rev(levels(data_long$variable)) 
    ) + 
    labs(y = "SHAP value (impact on model output)", x = "", color = "Feature value") 
  return(plot1)
}

# a function to standardize feature values into same range
std1 <- function(x){
  return ((x - min(x, na.rm = T))/(max(x, na.rm = T) - min(x, na.rm = T)))
}


# prep shap data
shap.prep <- function(shap  = shap_result, X_train = mydata$train_mm, top_n){
  require(ggforce)
  # descending order
  if (missing(top_n)) top_n <- dim(X_train)[2] # by default, use all features
  if (!top_n%in%c(1:dim(X_train)[2])) stop('supply correct top_n')
  require(data.table)
  shap_score_sub <- as.data.table(shap$shap_score)
  shap_score_sub <- shap_score_sub[, names(shap$mean_shap_score)[1:top_n], with = F]
  shap_score_long <- melt.data.table(shap_score_sub, measure.vars = colnames(shap_score_sub))
  
  # feature values: the values in the original dataset
  fv_sub <- as.data.table(X_train)[, names(shap$mean_shap_score)[1:top_n], with = F]
  # standardize feature values
  fv_sub_long <- melt.data.table(fv_sub, measure.vars = colnames(fv_sub))
  fv_sub_long[, stdfvalue := std1(value), by = "variable"]
  # SHAP value: value
  # raw feature value: rfvalue; 
  # standarized: stdfvalue
  names(fv_sub_long) <- c("variable", "rfvalue", "stdfvalue" )
  shap_long2 <- cbind(shap_score_long, fv_sub_long[,c('rfvalue','stdfvalue')])
  shap_long2[, mean_value := mean(abs(value)), by = variable]
  setkey(shap_long2, variable)
  return(shap_long2) 
}

var_importance <- function(shap_result, top_n = 10) {
  # Extract feature names and SHAP values
  feature_names <- names(shap_result$mean_shap_score)
  shap_values <- shap_result$mean_shap_score
  
  # Create a data frame with feature names and SHAP values
  shap_df <- tibble(var = feature_names, importance = shap_values)
  
  # Remove rows with NA values
  shap_df <- shap_df[complete.cases(shap_df), ]
  
  # Take the top_n features
  top_features <- shap_df %>%
    arrange(desc(abs(importance))) %>%
    head(top_n)
  
  return(top_features)
}

# Modified function to return top 10 SHAP values as a data frame
top_n_shap_df <- function(shap_result, top_n = 10) {
  feature_names <- names(shap_result$mean_shap_score)
  shap_values <- shap_result$mean_shap_score
  shap_df <- tibble(var = feature_names, importance = shap_values)
  
  # Remove rows with NA values
  shap_df <- shap_df[complete.cases(shap_df), ]
  
  # Take the top_n features
  top_features <- shap_df %>%
    arrange(desc(abs(importance))) %>%
    head(top_n)
  
  return(top_features)
}

#shap_values=predict(best.boost, dtrain_mQTL, predcontrib = TRUE, approxcontrib = F)

shap_result_bike = shap.score.rank(xgb_model = best.boost, 
                                   X_train =dtrain_mQTL,
                                   shap_approx = F
)

# Get the top 10 SHAP values as a data frame
top_10_shap_df <- top_n_shap_df(shap_result_bike, top_n = 10)
print(top_10_shap_df)

# Plot var importance based on SHAP
ggplot(top_10_shap_df, aes(x = reorder(var, importance), y = importance)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  theme_light() + 
  theme(axis.title.y = element_blank()) +
  labs(y = "Importance", x = "Feature")

##################################

###################################
#                                 #
#      Probe 5 - cg15597984       #
#                                 #            
###################################

# best model: Elastic net Ind GWAS before BH genetics + non-genetic
GenotypeFile <- read.csv("GWAS_trainset_beforeBH_indiv_probe5_snps_MAF10_genotype.csv",na = c("", "NA", "N/A",0))
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

# Define the columns to exclude
exclude_cols <- c("rs131969")

# Initialize a flag to indicate whether a row has been moved
row_moved <- FALSE

# Check each excluded column for "AA"
for (col in exclude_cols) {
  rows_to_move <- which(test.set[, col] == "A A")
  
  
  # If there are rows with "TT" or "CC" in the current column
  if (length(rows_to_move) > 0) {
    # Select a random row
    row_index <- sample(length(rows_to_move), 1)
    row_to_move <- rows_to_move[row_index]
    
    # Add the selected row to train.set
    train.set <- rbind(train.set, test.set[row_to_move, ])
    
    # Remove the selected row from test.set
    test.set <- test.set[-row_to_move, ]
    
    # Set the flag to indicate that a row has been moved
    row_moved <- TRUE
    
    # Break out of the loop to move only one row per iteration
    break
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
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 184))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 

#                     actual column + 1
new_row[, 3] <- "G G" #2 + 1
new_row[, 4] <- "A A" #3 + 1
new_row[, 142] <- "A A" #141 + 1
new_row[, 153] <- "A A" #152 + 1


colnames(new_row) <- colnames(train.set)

exclude_cols <- c(3,4,142,153)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:184, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
train.set <- rbind(train.set, new_row)

## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(seq_len(171), function(col_idx) length(unique(test.set[[col_idx]])) == 2))



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

write.csv(train.set,"Probe5_elasticnet_top_snps_trainset.csv", row.names = FALSE)


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- 0
lambda.vec <- 0

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

################# DALEX ########################

library(DALEX)

explainer <- explain(elastic, data = train.set[, -which(names(train.set) == "cg15597984")], y = train.set$cg15597984)

# Calculate feature importance using variable importance profiles (VIP)
vip_data <- variable_importance(explainer)

plot_data <- vip_data[, c("variable", "dropout_loss")]
plot_data <- plot_data[order(-plot_data$dropout_loss), ]

top_40_plot_data <- head(plot_data, 42)

ggplot(top_40_plot_data, aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Dropout Loss", title = "Top 10 Variable Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

############################# SHAP for elastic net ########################
library(ggplot2)

elastic_coef <- coef(elastic$finalModel, s = 0) # change s = 0 to your best lambda value

# Assuming 'elastic_coef' is the coefficient matrix obtained from the Elastic Net model
feature_importance <- data.frame(
  Feature = rownames(elastic_coef),
  Coefficient = as.numeric(elastic_coef)
)

# Take absolute values of the coefficients for importance
feature_importance$AbsoluteCoefficient <- abs(feature_importance$Coefficient)

# Remove the intercept (if it exists)
feature_importance <- feature_importance[feature_importance$Feature != "(Intercept)",]

# Sort the features by their absolute coefficients in descending order
feature_importance <- feature_importance[order(-feature_importance$AbsoluteCoefficient),]

# Select the top 10 features
top_10_features <- head(feature_importance, 10)

# Create a bar plot for the top 10 features
ggplot(top_10_features, aes(x = reorder(Feature, AbsoluteCoefficient), y = AbsoluteCoefficient)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Feature", y = "Importance") +
  theme_minimal()



###################################
#                                 #
#      Probe 7 - cg20046859       #
#                                 #            
###################################

# best model: Elastic net Ind GWAS uncorrected genetics only

GenotypeFile <- read.csv("GWAS_trainset_beforeBH_indiv_probe7_snps_MAF10_genotype.csv",na = c("", "NA", "N/A",0))
Probe <- read.csv("Probe7_cg20046859.csv")
standard.set <- read.csv("standard_testset_FID_v3.csv")

#merge the files by FID
MainData <- merge(GenotypeFile, Probe, by = "FID")

test.set <- merge(MainData, standard.set, by = "FID")
train.set <- MainData[!(MainData$FID %in% test.set$FID), ]

#remove ID, sex and phenotype columns
train.set <- train.set[, -c(1:6)]
test.set <- test.set[, -c(1:6)]

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

# for if test.set has 1 more genotype than train.set
# select row with that extra genotype in that column -> move it to train.set & remove from test.set

# Compare the number of categories in each column of train.set and test.set
for (col_name in colnames(train.set)) {
  if (length(unique(train.set[[col_name]])) < length(unique(test.set[[col_name]]))) {
    diff_categories <- setdiff(unique(test.set[[col_name]]), unique(train.set[[col_name]]))
    cat(paste0("Column '", col_name, "' has more categories in df2:", paste(diff_categories, collapse = ", "), "\n"))
  }
}


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
new_row <- data.frame(matrix(NA, nrow = 1, ncol = 154))


## to check: table(train.set$X.82)
## the letters are the ones with only 1 category in that col 

#                     actual column + 1
new_row[, 122] <- "G G" #121 + 1
new_row[, 143] <- "A A" #142 + 1



colnames(new_row) <- colnames(train.set)

exclude_cols <- c(122,143)
# Loop over each column, except the specific one, and randomly assign a value based on the existing classes
for (col in setdiff(1:154, exclude_cols)) {
  existing_classes <- unique(train.set[, col])
  new_row[, col] <- existing_classes[ceiling(runif(1, min = 1, max = length(existing_classes)))]
}

# Add the new row to the original dataframe
train.set <- rbind(train.set, new_row)


## step 3
# Find columns in test.set that have only 2 unique values
cols_with_2_classes <- which(sapply(test.set, function(x) length(unique(x))) == 2)

# Loop through the columns with 2 classes
new_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(test.set)))

for (col in cols_with_2_classes) {
  print(col)
  # Find the missing class in df1
  missing_class <- setdiff(levels(factor(train.set[[col]])), levels(factor(test.set[[col]])))
  
  
  #colnames(new_row) <- colnames(test.set)
  new_row[[col]] <- missing_class
}

# Add random values to the rest of the columns in the new row
for (i in 1:ncol(test.set)) {
  if (!(i %in% cols_with_2_classes)) {
    unique_vals <- unique(test.set[[i]])
    new_val <- unique_vals[ceiling(runif(1, min=1, max=length(unique_vals)))]
    new_row[[i]] <- new_val
  }
}

# Add the new row to test.set
colnames(new_row) <- colnames(test.set)
test.set <- rbind(test.set, new_row)

# Make sure all columns in test.set have 3 unique values
for (col in names(test.set)) {
  if (length(unique(test.set[[col]])) < 3) {
    stop("Not all columns in test.set have 3 unique values.")
  }
}

write.csv(train.set,"Probe7_elasticnet_top_snps_trainset.csv", row.names = FALSE)


#Training the Elastic Net model:

library(caret)
library(glmnet)

alpha.vec <- 0
lambda.vec <- 0

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

################# DALEX ########################

library(DALEX)

explainer <- explain(elastic, data = train.set[, -which(names(train.set) == "cg20046859")], y = train.set$cg20046859)

# Calculate feature importance using variable importance profiles (VIP)
vip_data <- variable_importance(explainer)

plot_data <- vip_data[, c("variable", "dropout_loss")]
plot_data <- plot_data[order(-plot_data$dropout_loss), ]

top_40_plot_data <- head(plot_data, 30)

ggplot(top_40_plot_data, aes(x = reorder(variable, dropout_loss), y = dropout_loss)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Dropout Loss", title = "Top 10 Variable Importance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()

############################# SHAP for elastic net ########################
library(ggplot2)

elastic_coef <- coef(elastic$finalModel, s = 0) # change s = 0 to your best lambda value

# Assuming 'elastic_coef' is the coefficient matrix obtained from the Elastic Net model
feature_importance <- data.frame(
  Feature = rownames(elastic_coef),
  Coefficient = as.numeric(elastic_coef)
)

# Take absolute values of the coefficients for importance
feature_importance$AbsoluteCoefficient <- abs(feature_importance$Coefficient)

# Remove the intercept (if it exists)
feature_importance <- feature_importance[feature_importance$Feature != "(Intercept)",]

# Sort the features by their absolute coefficients in descending order
feature_importance <- feature_importance[order(-feature_importance$AbsoluteCoefficient),]

# Select the top 10 features
top_10_features <- head(feature_importance, 10)

# Create a bar plot for the top 10 features
ggplot(top_10_features, aes(x = reorder(Feature, AbsoluteCoefficient), y = AbsoluteCoefficient)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(x = "Feature", y = "Importance") +
  theme_minimal()

