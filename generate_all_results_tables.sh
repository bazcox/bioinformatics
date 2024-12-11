#!/bin/bash
# generate_all_results_tables.sh
# Generates results tables for all combinations of feature sets, scenarios, and models.

echo "Generating tables for all combinations..."

feature_sets=("GWAS_beforeBH_combined" "GWAS_afterBH_combined" "GTEx")
scenarios=("genetic_only" "genetic_pcs_nongenetic")
models=("LinearRegression" "ElasticNet" "XGBoost")

# Loop through each combination and generate tables
for feature_set in "${feature_sets[@]}"; do
    for scenario in "${scenarios[@]}"; do
        for model in "${models[@]}"; do
            echo "Creating table for Feature_Set=${feature_set}, Scenario=${scenario}, Model=${model}"
            python generate_results_tables.py --feature_set "$feature_set" --scenario "$scenario" --model "$model"
        done
    done
done

echo "All tables generated."
