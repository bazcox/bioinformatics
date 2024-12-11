# CYP2D6 Methylation Prediction Replication

This project attempts to replicate and explore the methods from a research paper comparing different feature selection and machine learning approaches for predicting CYP2D6 methylation levels from genetic variation data. The original paper employed a variety of techniques including Linear Regression, Elastic Net, and XGBoost, as well as different feature sets (GWAS before/after BH correction, GTEx eQTLs, etc.), and considered both genetic features alone as well as additional principal components and non-genetic demographic variables.

Our replication involves:

* Data loading and preprocessing (genotypes, methylation, demographic, and PCA files)
* Splitting into training and test sets based on a standard test ID file
* Running multiple scenarios (e.g., genetic_only, genetic_pcs_nongenetic)
* Applying Linear Regression, Elastic Net, and XGBoost models
* Summarizing results in a manner similar to the paper (RMSE, R² for train/test)
* Generating tables in CSV and PNG image formats for easier comparison

## Project Overview

1. **run_model_pipeline.py**: Runs a single scenario (given a probe, feature set, and scenario) and appends results to a results_all.csv file. It trains Linear Regression, Elastic Net, and XGBoost, then prints and saves performance metrics (R² and MSE for train/test sets).

2. **execute_all_runs.sh**: A bash script that iterates over multiple probes, feature sets, and scenarios, calling run_model_pipeline.py for each combination. It produces results_all.csv containing all aggregated results.

3. **generate_all_results_tables.sh**: Given results_all.csv and arguments specifying a particular feature set, scenario, and model, it aggregates the results for the selected model and produces a table similar to those in the original paper. The table is printed to the console and also saved as a .png image.

4. **generate_all_tables.sh**: Runs generate_all_results_tables.sh for all combinations of feature sets, scenarios, and models, generating a series of tables in the tables directory.

## Setting Up the Environment

We recommend using [conda](https://anaconda.org/anaconda/conda) for managing environments:

```bash
# Create a new environment
conda create -n cyp2d6_env python=3.10 -y

# Activate the environment
conda activate cyp2d6_env

# Install necessary libraries
conda install pandas numpy scikit-learn matplotlib xgboost -y

# If needed, shap can be installed (optional):
conda install shap -y
```
Ensure that the data files are placed in the data/ directory and are organized as follows:

```scss
data/
 ├─ GUSTO_nongenetic_demo_withNA.csv
 ├─ genotypes/
 │   ├─ standard_testset_FID_v3.csv
 │   ├─ GTEx08snps_MAF10_genotype.csv
 │   ├─ GWAS_trainset_beforeBH_ALL_v2_genotype.csv
 │   ├─ GWAS_trainset_afterBH_ALL_v2_genotype.xlsx
 │   ... (other genotype files)
 ├─ methylation/
 │   ├─ Probe1_cg04692870_v2.csv
 │   ├─ Probe3_cg09322432.csv
 │   ├─ Probe5_cg15597984.csv
 │   ... (other probe files)
 ├─ pca/
 │   ├─ probe1_v2_genotype_subset_LDP_pca.eigenvec_5PCs_allsamples.csv
 │   ├─ probe3_genotype_subset_LDP_pca.eigenvec_5PCs.csv
 │   ├─ probe5_genotype_subset_LDP_pca.eigenvec_5PCs.csv
 │   ... (other PCA files)
```

## Example Workflow
1. **Generate Results**: After setting up your environment and placing all data in the data/ directory:

```bash
./execute_all_runs.sh
```
This will run run_model_pipeline.py over all combinations of probes, feature sets, and scenarios, producing results_all.csv.

2. **Generate Tables**: Once you have results_all.csv, you can create tables:
```bash
./generate_all_tables.sh
```
This script will call generate_all_results_tables.sh for each combination and store the resulting .png tables in ./tables.

## Running the Scripts

* **run_model_pipeline.py**:

```python 
python run_model_pipeline.py --probe Probe1 --feature_set GWAS_beforeBH_combined --scenario genetic_only
``` 
This will run the scenario for Probe1 with GWAS_beforeBH_combined features in the genetic_only scenario. It appends results to results_all.csv.

* **execute_all_runs.sh**:

```bash
./execute_all_runs.sh
```
Runs run_model_pipeline.py over all predefined combinations of probes, feature sets, and scenarios, generating results_all.csv.

* **generate_all_results_tables.sh**:

```python
python generate_all_results_tables.sh --feature_set GWAS_beforeBH_combined --scenario genetic_only --model LinearRegression
```
This reads from results_all.csv and generates a table (both console and PNG file) summarizing performance for the specified parameters.

* **generate_all_tables.sh**:

```bash
./generate_all_tables.sh
```
Iterates over all feature sets, scenarios, and models, producing multiple tables in ./tables.

## Further Adjustments

* To add or remove probes, edit the arrays in execute_all_runs.sh.
* To change models or parameters, modify run_model_pipeline.py or generate_all_results_tables.sh.
* If you need more complex hyperparameter tuning, adjust the parameter grids in run_model_pipeline.py for Elastic Net or XGBoost.