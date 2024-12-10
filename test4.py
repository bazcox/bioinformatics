import pandas as pd
import numpy as np
import os
import argparse
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV, LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV
import xgboost as xgb

# Attempt to import shap for feature importance analysis
# do_shap_analysis = True
do_shap_analysis = False
try:
    import shap
except ImportError:
    do_shap_analysis = False

parser = argparse.ArgumentParser(description="Run single scenario on specified probe and feature set.")
parser.add_argument("--probe", type=str, required=True)
parser.add_argument("--feature_set", type=str, required=True)
parser.add_argument("--scenario", type=str, required=True)
parser.add_argument("--output", type=str, default="results_all.csv", help="Output CSV file to store results.")
args = parser.parse_args()

probe_name = args.probe
feature_set_name = args.feature_set
scenario_name = args.scenario
output_file = args.output

data_dir = 'data'
geno_dir = os.path.join(data_dir, 'genotypes')
meth_dir = os.path.join(data_dir, 'methylation')
pca_dir = os.path.join(data_dir, 'pca')

probes = {
    'Probe1': ('Probe1_cg04692870_v2.csv', 'cg04692870'),
    'Probe2': ('Probe2_cg07016288.csv', 'cg07016288'),
    'Probe3': ('Probe3_cg09322432.csv', 'cg09322432'),
    'Probe4': ('Probe4_cg10840135.csv', 'cg10840135'),
    'Probe5': ('Probe5_cg15597984.csv', 'cg15597984'),
    'Probe6': ('Probe6_cg17498424_v2.csv', 'cg17498424'),
    'Probe7': ('Probe7_cg20046859.csv', 'cg20046859'),
    'Probe8': ('Probe8_cg22650942_v2.csv', 'cg22650942')
}

pca_files = {
    'Probe1': 'probe1_v2_genotype_subset_LDP_pca.eigenvec_5PCs_allsamples.csv',
    'Probe3': 'probe3_genotype_subset_LDP_pca.eigenvec_5PCs.csv',
    'Probe4': 'probe4_genotype_subset_LDP_pca.eigenvec_5PCs.csv',
    'Probe5': 'probe5_genotype_subset_LDP_pca.eigenvec_5PCs.csv',
    'Probe7': 'probe7_genotype_subset_LDP_pca.eigenvec_5PCs.csv'
}

feature_sets = {
    'GWAS_beforeBH_combined': {
        'description': 'GWAS before BH correction (combined)',
        'genotype_file': 'GWAS_trainset_beforeBH_ALL_v2_genotype.csv'
    },
    'GWAS_afterBH_combined': {
        'description': 'GWAS after BH correction (combined)',
        'genotype_file': 'GWAS_trainset_afterBH_ALL_v2_genotype.xlsx'
    },
    'GTEx': {
        'description': 'GTEx eQTLs',
        'genotype_file': 'GTEx08snps_MAF10_genotype.csv'
    }
}

scenarios = {
    'genetic_only': {'use_demographics': False, 'use_pcs': False},
    'genetic_pcs_nongenetic': {'use_demographics': True, 'use_pcs': True}
}

if probe_name not in probes:
    raise ValueError(f"Probe {probe_name} not found.")
if feature_set_name not in feature_sets:
    raise ValueError(f"Feature set {feature_set_name} not found.")
if scenario_name not in scenarios:
    raise ValueError(f"Scenario {scenario_name} not found.")

scenario = scenarios[scenario_name]
(use_demographics, use_pcs) = (scenario['use_demographics'], scenario['use_pcs'])

demo_df = pd.read_csv(os.path.join(data_dir, 'GUSTO_nongenetic_demo_withNA.csv'))
demo_df = demo_df.rename(columns={'masked ids': 'sample_id'})
demo_df = demo_df.dropna(subset=['sample_id']).drop_duplicates('sample_id')

demo_features = [
    'mother_income_cat', 'household_income_cat', 'accommodation_cat',
    'Child_ethnicity_merged', 'mother_highest_education_cat',
    'child.s_sex_binary', 'mother_age_recruitment'
]

test_ids_df = pd.read_csv(os.path.join(geno_dir, 'standard_testset_FID_v3.csv'))
test_ids_df = test_ids_df.rename(columns={'Masked IDs': 'sample_id'})
test_ids = test_ids_df['sample_id'].unique()
test_ids_set = set(test_ids)

all_ids = demo_df['sample_id'].unique()
train_ids_set = set(all_ids) - test_ids_set

meth_file, meth_col = probes[probe_name]
meth_df = pd.read_csv(os.path.join(meth_dir, meth_file))
id_col = [c for c in meth_df.columns if 'masked' in c.lower()]
if len(id_col) != 1:
    raise ValueError(f"Unique masked column not found in {meth_file}")
meth_df = meth_df.rename(columns={id_col[0]: 'sample_id'})
meth_df = meth_df.dropna(subset=['sample_id']).drop_duplicates('sample_id')

merged_df = demo_df.merge(meth_df[['sample_id', meth_col]], on='sample_id', how='inner')

if probe_name in pca_files:
    pca_df = pd.read_csv(os.path.join(pca_dir, pca_files[probe_name]))
    if 'Masked IDs' in pca_df.columns:
        pca_df = pca_df.rename(columns={'Masked IDs': 'sample_id'})
    else:
        raise ValueError("No 'Masked IDs' column found in PCA file.")
    pcs = [c for c in pca_df.columns if 'PC' in c.upper()]
    if 'sample_id' not in pca_df.columns:
        raise ValueError("No 'sample_id' column found after renaming in PCA file.")
    pca_df = pca_df[['sample_id'] + pcs]
    merged_df = merged_df.merge(pca_df, on='sample_id', how='left')
else:
    pcs = []

train_df = merged_df[merged_df['sample_id'].isin(train_ids_set)].copy()
test_df = merged_df[merged_df['sample_id'].isin(test_ids_set)].copy()

geno_file = feature_sets[feature_set_name]['genotype_file']
if geno_file.endswith('.xlsx'):
    geno_df = pd.read_excel(os.path.join(geno_dir, geno_file))
else:
    geno_df = pd.read_csv(os.path.join(geno_dir, geno_file), low_memory=False)

if 'masked' in geno_df.columns:
    geno_df = geno_df.rename(columns={'masked': 'sample_id'})

geno_df = geno_df.dropna(subset=['sample_id']).drop_duplicates('sample_id')

exclude_cols = ['sample_id', 'PID', 'MID', 'Sex', 'Phenotype']
snp_cols = [c for c in geno_df.columns if c not in exclude_cols]

train_full = train_df.merge(geno_df[['sample_id']+snp_cols], on='sample_id', how='inner')
test_full = test_df.merge(geno_df[['sample_id']+snp_cols], on='sample_id', how='inner')

base_features = snp_cols[:]

if use_demographics:
    for col in demo_features:
        if col not in train_full.columns:
            print(f"Warning: {col} not found in merged data.")
        else:
            base_features.append(col)

if use_pcs and pcs:
    base_features.extend(pcs)

base_features = list(dict.fromkeys(base_features))

y_train = train_full[meth_col]
X_train = train_full[base_features]
y_test = test_full[meth_col]
X_test = test_full[base_features]

cat_cols = [c for c in X_train.columns if X_train[c].dtype == 'object']
if cat_cols:
    X_train = pd.get_dummies(X_train, columns=cat_cols, drop_first=True)
    X_test = pd.get_dummies(X_test, columns=cat_cols, drop_first=True)
    X_test = X_test.reindex(columns=X_train.columns, fill_value=0)

imputer = SimpleImputer(strategy='mean')
X_train_imputed = imputer.fit_transform(X_train)
X_test_imputed = imputer.transform(X_test)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_imputed)
X_test_scaled = scaler.transform(X_test_imputed)

# LinearRegression (no regularization)
lr = LinearRegression()
lr.fit(X_train_scaled, y_train)
y_pred_lr_train = lr.predict(X_train_scaled)
y_pred_lr_test = lr.predict(X_test_scaled)
train_mse_lr = mean_squared_error(y_train, y_pred_lr_train)
train_r2_lr = r2_score(y_train, y_pred_lr_train)
mse_lr = mean_squared_error(y_test, y_pred_lr_test)
r2_lr = r2_score(y_test, y_pred_lr_test)

# Elastic Net
en = ElasticNetCV(l1_ratio=[0.1,0.5,0.9],
                  alphas=[0.001,0.01,0.1,1,10],
                  cv=3,
                  max_iter=10000)
en.fit(X_train_scaled, y_train)
y_pred_en_train = en.predict(X_train_scaled)
y_pred_en_test = en.predict(X_test_scaled)
train_mse_en = mean_squared_error(y_train, y_pred_en_train)
train_r2_en = r2_score(y_train, y_pred_en_train)
mse_en = mean_squared_error(y_test, y_pred_en_test)
r2_en = r2_score(y_test, y_pred_en_test)

# XGBoost
param_grid = {
    'n_estimators':[100],
    'max_depth':[3],
    'learning_rate':[0.1]
}
xg_model = xgb.XGBRegressor(random_state=42, eval_metric='rmse')
gs = GridSearchCV(xg_model, param_grid, cv=3, scoring='neg_mean_squared_error')
gs.fit(X_train_scaled, y_train)
best_xg = gs.best_estimator_
y_pred_xg_train = best_xg.predict(X_train_scaled)
y_pred_xg_test = best_xg.predict(X_test_scaled)
train_mse_xg = mean_squared_error(y_train, y_pred_xg_train)
train_r2_xg = r2_score(y_train, y_pred_xg_train)
mse_xg = mean_squared_error(y_test, y_pred_xg_test)
r2_xg = r2_score(y_test, y_pred_xg_test)

print("\n========== RESULTS ==========")
print(f"Probe: {probe_name}")
print(f"Feature Set: {feature_set_name}")
print(f"Scenario: {scenario_name}")

print("\nLinear Regression Performance:")
print(f"  R² train: {train_r2_lr:.4f}")
print(f"  MSE train: {train_mse_lr:.4f}")
print(f"  R² test: {r2_lr:.4f}")
print(f"  MSE test: {mse_lr:.4f}")

print("\nElastic Net Performance:")
print(f"  R² train: {train_r2_en:.4f}")
print(f"  MSE train: {train_mse_en:.4f}")
print(f"  R² test: {r2_en:.4f}")
print(f"  MSE test: {mse_en:.4f}")

print("\nXGBoost Performance:")
print(f"  R² train: {train_r2_xg:.4f}")
print(f"  MSE train: {train_mse_xg:.4f}")
print(f"  R² test: {r2_xg:.4f}")
print(f"  MSE test: {mse_xg:.4f}")

if do_shap_analysis:
    explainer = shap.TreeExplainer(best_xg)
    shap_values = explainer.shap_values(X_test_scaled)
    feature_names = X_test.columns
    mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
    top5_idx = np.argsort(mean_abs_shap)[-5:]
    top_features = feature_names[top5_idx]

    print("\nTop 5 Features by Mean Absolute SHAP Value:")
    for f in reversed(top_features):
        print(f"  {f}")
    shap.summary_plot(shap_values, X_test_scaled, feature_names=feature_names)

# Append results to CSV
results = {
    'Probe': probe_name,
    'Feature_Set': feature_set_name,
    'Scenario': scenario_name,
    'LinearRegression_R2_train': train_r2_lr,
    'LinearRegression_MSE_train': train_mse_lr,
    'LinearRegression_R2_test': r2_lr,
    'LinearRegression_MSE_test': mse_lr,
    'ElasticNet_R2_train': train_r2_en,
    'ElasticNet_MSE_train': train_mse_en,
    'ElasticNet_R2_test': r2_en,
    'ElasticNet_MSE_test': mse_en,
    'XGBoost_R2_train': train_r2_xg,
    'XGBoost_MSE_train': train_mse_xg,
    'XGBoost_R2_test': r2_xg,
    'XGBoost_MSE_test': mse_xg
}

df_out = pd.DataFrame([results])
if not os.path.exists(output_file):
    df_out.to_csv(output_file, index=False)
else:
    df_out.to_csv(output_file, mode='a', index=False, header=False)
