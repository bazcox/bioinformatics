import pandas as pd
import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.decomposition import PCA
from sklearn.feature_selection import RFE
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from xgboost import XGBRegressor
import shap
import os

# ---------------------------
# Load and Clean Test Dataset
# ---------------------------
test_file = os.path.join(
    "src", "Updated linear regression files and R codes",
    "GWAS_75testset_all_afterBH_probe5_v2_converted.csv"
)

test_df = pd.read_csv(test_file).dropna()

# ---------------------------
# Drop Non-Numeric Columns
# ---------------------------
non_numeric_cols = test_df.select_dtypes(include=['object']).columns
print(f"Non-Numeric Columns: {non_numeric_cols}")

# Drop non-numeric columns
test_df = test_df.drop(columns=non_numeric_cols, errors='ignore')

# ---------------------------
# Create New Features (Interaction Terms)
# ---------------------------
significant_snps = ['rs5758165_cat', 'rs114309992_cat', 'rs76550409_cat', 'rs5751117_cat', 'rs134906_cat']

for snp in significant_snps:
    for pc in ['PC1.x', 'PC2.x', 'PC3.x', 'PC4', 'PC5']:
        test_df[f'{snp}_{pc}_interaction'] = test_df[snp] * test_df[pc]

# Create dummy variables for categorical features
test_df['mother_chinesevmalay'] = np.where(test_df['mother_ethnicity_cat'].isin([1, 3]), 0, (test_df['mother_ethnicity_cat'] == 2).astype(int))

test_df['mother_chinesevindian'] = np.where(test_df['mother_ethnicity_cat'].isin([1, 2]), 0, (test_df['mother_ethnicity_cat'] == 3).astype(int))

# ---------------------------
# Define Target and Features
# ---------------------------
target = 'cg15597984'
X = test_df.drop(columns=[target])
y = test_df[target]

# ---------------------------
# Standardize Numeric Features
# ---------------------------
numeric_features = X.select_dtypes(include=['float64', 'int64']).columns
scaler = StandardScaler()
X[numeric_features] = scaler.fit_transform(X[numeric_features])

# ---------------------------
# Feature Selection using RFE
# ---------------------------
ridge_for_rfe = RidgeCV(alphas=np.logspace(-6, 6, 13), cv=5)
rfe = RFE(estimator=ridge_for_rfe, n_features_to_select=30)
X_rfe = rfe.fit_transform(X, y)
selected_features = X.columns[rfe.support_]

print("\nSelected Features by RFE:")
print(selected_features)

# ---------------------------
# Dimensionality Reduction with PCA
# ---------------------------
pca = PCA(n_components=20)
X_pca = pca.fit_transform(X_rfe)

# ---------------------------
# Ridge Regression with Polynomial Features
# ---------------------------
poly = PolynomialFeatures(degree=2, interaction_only=True, include_bias=False)
X_poly = poly.fit_transform(X_pca)

ridge = RidgeCV(alphas=np.logspace(-6, 6, 13), cv=5)
ridge.fit(X_poly, y)

# Ridge Regression Predictions
y_pred_ridge = ridge.predict(X_poly)
rmse_ridge = np.sqrt(mean_squared_error(y, y_pred_ridge))
r2_ridge = r2_score(y, y_pred_ridge)

print("\nRidge Regression Results")
print(f"RMSE: {rmse_ridge:.6f}")
print(f"R²: {r2_ridge:.6f}")

# ---------------------------
# XGBoost with Hyperparameter Tuning
# ---------------------------
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [3, 4, 5],
    'learning_rate': [0.01, 0.05, 0.1],
    'subsample': [0.7, 0.8, 1.0]
}

xgb_model = XGBRegressor(random_state=42)
grid_search = GridSearchCV(xgb_model, param_grid, cv=5, scoring='r2')
grid_search.fit(X_rfe, y)

best_xgb = grid_search.best_estimator_
y_pred_xgb = best_xgb.predict(X_rfe)
rmse_xgb = np.sqrt(mean_squared_error(y, y_pred_xgb))
r2_xgb = r2_score(y, y_pred_xgb)

print("\nXGBoost Results with Best Hyperparameters")
print(f"RMSE: {rmse_xgb:.6f}")
print(f"R²: {r2_xgb:.6f}")

# ---------------------------
# SHAP Analysis for XGBoost
# ---------------------------
explainer = shap.Explainer(best_xgb, X_rfe)
shap_values = explainer(X_rfe)

print("\nTop 10 Features by SHAP Values:")
shap.summary_plot(shap_values, X[selected_features], max_display=10)

# ---------------------------
# Ensemble of Ridge Regression and XGBoost
# ---------------------------
y_pred_ensemble = (y_pred_ridge + y_pred_xgb) / 2
rmse_ensemble = np.sqrt(mean_squared_error(y, y_pred_ensemble))
r2_ensemble = r2_score(y, y_pred_ensemble)

print("\nEnsemble Model Results")
print(f"RMSE: {rmse_ensemble:.6f}")
print(f"R²: {r2_ensemble:.6f}")

# ---------------------------
# Export Predictions and Actual Values
# ---------------------------
output_file = os.path.join("src", "Updated linear regression files and R codes", "predictions_output.csv")
results_df = pd.DataFrame({
    'Actual': y,
    'Ridge_Predicted': y_pred_ridge,
    'XGBoost_Predicted': y_pred_xgb,
    'Ensemble_Predicted': y_pred_ensemble
})
results_df.to_csv(output_file, index=False)

print(f"\nPredictions and actual values exported to: {output_file}")
