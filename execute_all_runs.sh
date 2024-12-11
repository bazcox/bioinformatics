#!/bin/bash
# execute_all_runs.sh
# This script runs the modeling pipeline for all combinations of probes, feature sets, and scenarios.

echo "Starting all runs..."

# Define combinations of probes, feature sets, and scenarios
probes=("Probe1" "Probe3" "Probe5")
feature_sets=("GWAS_beforeBH_combined" "GWAS_afterBH_combined" "GTEx")
scenarios=("genetic_only" "genetic_pcs_nongenetic")

# Loop through each combination and run the pipeline
for probe in "${probes[@]}"; do
    for feature_set in "${feature_sets[@]}"; do
        for scenario in "${scenarios[@]}"; do
            echo "Running: Probe=${probe}, Feature_Set=${feature_set}, Scenario=${scenario}"
            python run_model_pipeline.py --probe "$probe" --feature_set "$feature_set" --scenario "$scenario"
        done
    done
done

echo "All runs completed."
