import pandas as pd
import os
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.linear_model import LinearRegression, ElasticNetCV, LassoCV
from sklearn.metrics import mean_squared_error, r2_score
import xgboost as xgb

# ---------------------------
# Load and Clean Demographic Data
# ---------------------------
data_dir = 'data'
demo_df = pd.read_csv(os.path.join(data_dir, 'GUSTO_nongenetic_demo_withNA.csv'))
demo_df = demo_df.rename(columns={'masked ids': 'sample_id'}).dropna(subset=['sample_id']).drop_duplicates(subset='sample_id')

# ---------------------------
# Load and Clean Methylation Data
# ---------------------------
meth_file = os.path.join(data_dir, 'methylation', 'Probe1_cg04692870_v2.csv')
meth_df = pd.read_csv(meth_file)
meth_df = meth_df.rename(columns={'Probe 1 masked': 'sample_id'}).dropna(subset=['sample_id']).drop_duplicates(subset='sample_id')

# ---------------------------
# Load and Clean Genotype Data
# ---------------------------
geno_file = os.path.join(data_dir, 'genotypes', 'GTEx08snps_MAF10_genotype.csv')
geno_df = pd.read_csv(geno_file).rename(columns={'masked': 'sample_id'}).dropna(subset=['sample_id']).drop_duplicates(subset='sample_id')

# ---------------------------
# Load PCA Data
# ---------------------------
pca_file = os.path.join(data_dir, 'pca', 'probe1_v2_genotype_subset_LDP_pca.eigenvec_5PCs_allsamples.csv')
pca_df = pd.read_csv(pca_file)
pca_df = pca_df.rename(columns={'Masked IDs': 'sample_id'}).dropna(subset=['sample_id'])

# ---------------------------
# Merge DataFrames
# ---------------------------
combined_df = (
    demo_df
    .merge(geno_df, on='sample_id', how='inner')
    .merge(meth_df, on='sample_id', how='inner')
    .merge(pca_df, on='sample_id', how='inner')
)

print("\nCombined DataFrame head:")
print(combined_df.head())
print("Number of samples after merging:", combined_df.shape[0])

# ---------------------------
# Prepare Features and Target
# ---------------------------
y = combined_df['cg04692870']
X = combined_df.drop(['cg04692870', 'sample_id'], axis=1, errors='ignore')

# ---------------------------
# Imputation and Encoding
# ---------------------------
# Identify numeric and categorical columns
numeric_cols = X.select_dtypes(include=['number']).columns
categorical_cols = X.select_dtypes(include=['object', 'category']).columns

print("\nNumeric columns:", numeric_cols)
print("Categorical columns:", categorical_cols)

# Define imputers
numeric_imputer = SimpleImputer(strategy='mean')
categorical_imputer = SimpleImputer(strategy='most_frequent')

# Preprocess the data
preprocessor = ColumnTransformer(
    transformers=[
        ('num', numeric_imputer, numeric_cols),
        ('cat', categorical_imputer, categorical_cols)
    ]
)

# Fit and transform the data
X_imputed = preprocessor.fit_transform(X)

# One-hot encode categorical data
categorical_encoder = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
X_encoded = categorical_encoder.fit_transform(X[categorical_cols])

# Combine numeric and encoded categorical data
X_final = np.hstack((X_imputed[:, :len(numeric_cols)], X_encoded))

# ---------------------------
# Scaling
# ---------------------------
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_final)

# ---------------------------
# Train-Test Split
# ---------------------------
X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, test_size=0.25, random_state=42)

# ---------------------------
# Feature Selection Using Lasso
# ---------------------------
lasso = LassoCV(cv=5, max_iter=10000)
lasso.fit(X_train, y_train)

# Get the non-zero coefficient features
important_features = lasso.coef_ != 0
print(f"Number of important features selected by Lasso: {np.sum(important_features)}")

# Reduce X to only important features
X_train_selected = X_train[:, important_features]
X_test_selected = X_test[:, important_features]

# ---------------------------
# Model Training and Evaluation
# ---------------------------

# Linear Regression
lr = LinearRegression()
lr.fit(X_train_selected, y_train)
y_pred_lr = lr.predict(X_test_selected)
print("\nLinear Regression MSE:", mean_squared_error(y_test, y_pred_lr), "R^2:", r2_score(y_test, y_pred_lr))

# Elastic Net with Increased Iterations
en = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9], alphas=[0.0001, 0.001, 0.01, 0.1, 1, 10], cv=5, max_iter=10000)
en.fit(X_train_selected, y_train)
y_pred_en = en.predict(X_test_selected)
print("Elastic Net MSE:", mean_squared_error(y_test, y_pred_en), "R^2:", r2_score(y_test, y_pred_en))

# XGBoost with Grid Search
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [3, 5, 7],
    'learning_rate': [0.01, 0.05, 0.1],
    'subsample': [0.8, 1.0]
}
gs = GridSearchCV(xgb.XGBRegressor(random_state=42), param_grid, cv=5, scoring='r2')
gs.fit(X_train_selected, y_train)

best_xg = gs.best_estimator_
y_pred_xg = best_xg.predict(X_test_selected)
print("XGBoost Best Params:", gs.best_params_)
print("XGBoost MSE:", mean_squared_error(y_test, y_pred_xg), "R^2:", r2_score(y_test, y_pred_xg))
