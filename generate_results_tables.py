"""
generate_results_tables.py

Generates summary tables from the aggregated results CSV file. The tables
include performance metrics (RMSE, R²) for different models, feature sets, and scenarios.

Arguments:
    --feature_set: Feature set to analyze (e.g., 'GWAS_beforeBH_combined').
    --scenario: Scenario to analyze (e.g., 'genetic_only').
    --model: Model to analyze (e.g., 'LinearRegression').

Example:
    python generate_results_tables.py --feature_set GWAS_beforeBH_combined --scenario genetic_only --model XGBoost
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

# ---------------------------
# Parse Arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Create tables for a given Feature_Set and Scenario")
parser.add_argument("--feature_set", type=str, required=True, help="Feature set name")
parser.add_argument("--scenario", type=str, required=True, help="Scenario name")
parser.add_argument("--model", type=str, default="LinearRegression", help="Which model to use (LinearRegression, ElasticNet, XGBoost)")
parser.add_argument("--output_prefix", type=str, default="table", help="Prefix for output image file")
args = parser.parse_args()

chosen_feature_set = args.feature_set
chosen_scenario = args.scenario
model = args.model
output_prefix = args.output_prefix

# ---------------------------
# Load the Data
# ---------------------------
df = pd.read_csv("results_all.csv")

# Filter Data for the Chosen Feature Set and Scenario
subset = df[(df['Feature_Set'] == chosen_feature_set) & (df['Scenario'] == chosen_scenario)]

# Check if the subset is empty
if subset.empty:
    raise ValueError(f"No data found for Feature_Set={chosen_feature_set} and Scenario={chosen_scenario}")

# ---------------------------
# Define Column Names for the Metrics
# ---------------------------
probes = subset['Probe'].unique()
probes_list = sorted(probes)
final_columns = probes_list + ["Mean (SD)"]

# Construct column names based on the chosen model
r2_train_col = f"{model}_R2_train"
mse_train_col = f"{model}_MSE_train"
r2_test_col = f"{model}_R2_test"
mse_test_col = f"{model}_MSE_test"

# ---------------------------
# Helper Functions
# ---------------------------
def pivot_metric(metric_col):
    """Pivot the subset DataFrame for a given metric column."""
    # Create a copy of the subset to avoid SettingWithCopyWarning
    subset_copy = subset[['Scenario', 'Probe', metric_col]].copy()

    # Check for duplicates and remove them from the copy
    if subset_copy.duplicated().any():
        # print("Warning: Duplicate entries detected. Removing duplicates.")
        subset_copy.drop_duplicates(subset=['Scenario', 'Probe'], inplace=True)

    # Check for missing values and remove rows with NaN
    if subset_copy.isnull().any().any():
        # print(f"Warning: Missing values detected in {metric_col}. Filling with NaN.")
        subset_copy.dropna(subset=[metric_col], inplace=True)

    pivot = subset_copy.pivot(index='Scenario', columns='Probe', values=metric_col)
    pivot = pivot.iloc[0]
    return pivot

def mean_sd(series):
    """Compute mean and standard deviation and format as a string."""
    mean_val = series.mean()
    std_val = series.std()
    return f"{mean_val:.4f} ({std_val:.4f})"

def add_mean_sd_col(series):
    """Add mean and standard deviation column to the series."""
    vals = [f"{v:.4f}" for v in series]
    vals.append(mean_sd(series))
    return vals

# ---------------------------
# Compute Metrics
# ---------------------------
rmse_train = pivot_metric(mse_train_col).apply(np.sqrt)
r2_train = pivot_metric(r2_train_col)
rmse_test = pivot_metric(mse_test_col).apply(np.sqrt)
r2_test = pivot_metric(r2_test_col)

rmse_train_values = add_mean_sd_col(rmse_train)
r2_train_values = add_mean_sd_col(r2_train)
rmse_test_values = add_mean_sd_col(rmse_test)
r2_test_values = add_mean_sd_col(r2_test)

# ---------------------------
# Create Final Table
# ---------------------------
all_columns = ["Metric"] + final_columns
final_table = pd.DataFrame(index=range(4), columns=all_columns)

final_table.loc[0, "Metric"] = "RMSE train set"
final_table.loc[0, final_columns] = rmse_train_values

final_table.loc[1, "Metric"] = "R² train set"
final_table.loc[1, final_columns] = r2_train_values

final_table.loc[2, "Metric"] = "RMSE test set"
final_table.loc[2, final_columns] = rmse_test_values

final_table.loc[3, "Metric"] = "R² test set"
final_table.loc[3, final_columns] = r2_test_values

print(final_table.to_string(index=False))

# ---------------------------
# Save Table as Image with Header
# ---------------------------
fig, ax = plt.subplots(figsize=(10, 4))  # Increased height to accommodate the header
ax.axis('tight')
ax.axis('off')

# Add the header text above the table
header_text = f"Feature Set: {chosen_feature_set} | Scenario: {chosen_scenario} | Model: {model}"
fig.text(0.5, 0.95, header_text, fontsize=14, ha='center', va='center', fontweight='bold')

# Create the table
table_data = [final_table.columns.tolist()] + final_table.values.tolist()
the_table = ax.table(cellText=table_data,
                     colLabels=None,
                     loc='center',
                     cellLoc='center')
the_table.set_fontsize(10)
the_table.scale(1, 2)

# Adjust layout to leave space for the header
plt.subplots_adjust(top=0.85, bottom=0.05)

# Ensure output directory exists
output_dir = "tables"
os.makedirs(output_dir, exist_ok=True)

image_filename = f"{output_prefix}_{chosen_feature_set}_{chosen_scenario}_{model}.png"
image_path = os.path.join(output_dir, image_filename)
plt.savefig(image_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"Table image saved as {image_filename}\n")


