import pandas as pd
import numpy as np
import os

from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV
import xgboost as xgb

# ===========================================
# Paths
# ===========================================
data_dir = 'data'
geno_dir = os.path.join(data_dir, 'genotypes')
meth_dir = os.path.join(data_dir, 'methylation')
pca_dir = os.path.join(data_dir, 'pca')

# ===========================================
# Probes and Files
# ===========================================
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

# ===========================================
# Feature Sets (without the large 2MB_range set)
# ===========================================
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

# ===========================================
# Load Demographics
# ===========================================
demo_df = pd.read_csv(os.path.join(data_dir, 'GUSTO_nongenetic_demo_withNA.csv'))
demo_df = demo_df.rename(columns={'masked ids': 'sample_id'})
demo_df = demo_df.dropna(subset=['sample_id']).drop_duplicates('sample_id')

demo_features = [
    'mother_income_cat', 'household_income_cat', 'accommodation_cat',
    'Child_ethnicity_merged', 'mother_highest_education_cat',
    'child.s_sex_binary', 'mother_age_recruitment'
]

# ===========================================
# Standard Test Set
# ===========================================
test_ids_df = pd.read_csv(os.path.join(geno_dir, 'standard_testset_FID_v3.csv'))
test_ids_df = test_ids_df.rename(columns={'Masked IDs': 'sample_id'})
test_ids = test_ids_df['sample_id'].unique()
test_ids_set = set(test_ids)

all_ids = demo_df['sample_id'].unique()
train_ids_set = set(all_ids) - test_ids_set

# ===========================================
# Scenarios
# ===========================================
scenarios = [
    {'name': 'genetic_only', 'use_demographics': False, 'use_pcs': False},
    {'name': 'genetic_pcs_nongenetic', 'use_demographics': True, 'use_pcs': True}
]

results = []

for probe_name, (meth_file, meth_col) in probes.items():
    meth_df = pd.read_csv(os.path.join(meth_dir, meth_file))
    id_col = [c for c in meth_df.columns if 'masked' in c.lower()]
    if len(id_col) != 1:
        print(f"Warning: Could not identify unique masked column in {meth_file}")
        continue
    meth_df = meth_df.rename(columns={id_col[0]: 'sample_id'})
    meth_df = meth_df.dropna(subset=['sample_id']).drop_duplicates('sample_id')
    
    merged_df = demo_df.merge(meth_df[['sample_id', meth_col]], on='sample_id', how='inner')
    
    if probe_name in pca_files:
        pca_df = pd.read_csv(os.path.join(pca_dir, pca_files[probe_name]))
        if 'Masked IDs' in pca_df.columns:
            pca_df = pca_df.rename(columns={'Masked IDs': 'sample_id'})
        else:
            print(f"PCA file {pca_files[probe_name]} columns:", pca_df.columns)
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
    
    for fs_name, fs_info in feature_sets.items():
        geno_file = fs_info['genotype_file']
        
        if geno_file.endswith('.xlsx'):
            geno_df = pd.read_excel(os.path.join(geno_dir, geno_file))
        else:
            # Use low_memory=False to reduce dtype warnings and speed parsing
            geno_df = pd.read_csv(os.path.join(geno_dir, geno_file), low_memory=False)
        
        if 'masked' in geno_df.columns:
            geno_df = geno_df.rename(columns={'masked': 'sample_id'})
        
        geno_df = geno_df.dropna(subset=['sample_id']).drop_duplicates('sample_id')
        
        exclude_cols = ['sample_id', 'PID', 'MID', 'Sex', 'Phenotype']
        snp_cols = [c for c in geno_df.columns if c not in exclude_cols]
        
        train_full = train_df.merge(geno_df[['sample_id']+snp_cols], on='sample_id', how='inner')
        test_full = test_df.merge(geno_df[['sample_id']+snp_cols], on='sample_id', how='inner')
        
        for scenario in scenarios:
            scenario_name = scenario['name']
            
            base_features = snp_cols[:]
            
            if scenario['use_demographics']:
                for col in demo_features:
                    if col not in train_full.columns:
                        print(f"Warning: {col} not found in merged data.")
                    else:
                        base_features.append(col)
            
            if scenario['use_pcs'] and pcs:
                base_features.extend(pcs)
            
            base_features = list(dict.fromkeys(base_features))
            
            y_train = train_full[meth_col]
            X_train = train_full[base_features]
            y_test = test_full[meth_col]
            X_test = test_full[base_features]
            
            # If SNPs or demographics are categorical, encode
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
            
            # Increase max_iter and reduce cv in ElasticNetCV
            # Increase max_iter and set a stricter tolerance
            en = ElasticNetCV(l1_ratio=[0.1, 0.5, 0.9],
                            alphas=[0.001,0.01,0.1,1,10],
                            cv=3,
                            max_iter=50000,
                            tol=1e-6)
            en.fit(X_train_scaled, y_train)
            y_pred_en = en.predict(X_test_scaled)
            mse_en = mean_squared_error(y_test, y_pred_en)
            r2_en = r2_score(y_test, y_pred_en)
            
            # Reduce parameter grid for XGBoost for faster run
            param_grid = {
                'n_estimators':[100],
                'max_depth':[3],
                'learning_rate':[0.1]
            }
            xg_model = xgb.XGBRegressor(random_state=42, eval_metric='rmse')
            gs = GridSearchCV(xg_model, param_grid, cv=3, scoring='neg_mean_squared_error')
            gs.fit(X_train_scaled, y_train)
            best_xg = gs.best_estimator_
            y_pred_xg = best_xg.predict(X_test_scaled)
            mse_xg = mean_squared_error(y_test, y_pred_xg)
            r2_xg = r2_score(y_test, y_pred_xg)
            
            results.append({
                'Probe': probe_name,
                'Feature_Set': fs_name,
                'Scenario': scenario_name,
                'ElasticNet_R2': r2_en,
                'ElasticNet_MSE': mse_en,
                'XGBoost_R2': r2_xg,
                'XGBoost_MSE': mse_xg
            })

results_df = pd.DataFrame(results)
results_df.to_csv("model_results_comparison.csv", index=False)
print("Results saved to model_results_comparison.csv")
print(results_df.head())
