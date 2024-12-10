import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the aggregated results CSV
df = pd.read_csv("results_all.csv")

# Example: replicate a table similar to the paper
# Let's say we want to get the mean and std of test R² for ElasticNet across probes, grouped by Feature_Set.
# Modify these groupings as needed:
group_cols = ['Probe', 'Feature_Set', 'Scenario']

# Aggregate mean and std for each model's test R² and MSE
agg_metrics = {
    'LinearRegression_R2_test': ['mean','std'],
    'LinearRegression_MSE_test': ['mean','std'],
    'ElasticNet_R2_test': ['mean','std'],
    'ElasticNet_MSE_test': ['mean','std'],
    'XGBoost_R2_test': ['mean','std'],
    'XGBoost_MSE_test': ['mean','std']
}

summary = df.groupby(group_cols).agg(agg_metrics)

# Print the summary table
print("\n=== Aggregated Test Performance ===")
print(summary)

# If you want a table similar to the paper that averages over Probes, you can drop 'Probe' from the group:
# For example:
summary_by_feature_scenario = df.groupby(['Feature_Set','Scenario']).agg(agg_metrics)
print("\n=== Mean and Std by Feature_Set and Scenario ===")
print(summary_by_feature_scenario)

# -------------------------
# Create a Graph Example
# -------------------------
# Suppose we want to plot the average test R² for ElasticNet across different Feature_Sets for a given scenario.
scenario_to_plot = 'genetic_only'  # choose a scenario
model_to_plot = 'ElasticNet_R2_test'

# Filter data for chosen scenario
subset = df[df['Scenario']==scenario_to_plot]

# Compute mean R² by Feature_Set (and maybe across all Probes)
mean_values = subset.groupby('Feature_Set')[model_to_plot].mean()

plt.figure(figsize=(8,6))
mean_values.plot(kind='bar', color='skyblue')
plt.title(f"Mean {model_to_plot} by Feature_Set ({scenario_to_plot})")
plt.xlabel("Feature_Set")
plt.ylabel(f"{model_to_plot}")
plt.axhline(y=0, color='black', linewidth=1)
plt.tight_layout()
plt.savefig("mean_R2_barplot.png")
plt.show()

# You can similarly create plots for other models or metrics as needed.
