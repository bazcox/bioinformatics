#!/usr/bin/env bash

PROBES=("Probe1" "Probe3" "Probe5")
FEATURE_SETS=("GWAS_beforeBH_combined" "GWAS_afterBH_combined" "GTEx")
SCENARIOS=("genetic_only" "genetic_pcs_nongenetic")

SCRIPT_DIR=$(dirname "$0")

# Remove old results file if exists
rm -f results_all.csv

echo "Starting all runs..."
for probe in "${PROBES[@]}"; do
    for feature_set in "${FEATURE_SETS[@]}"; do
        for scenario in "${SCENARIOS[@]}"; do
            echo "Running: Probe=$probe, Feature_Set=$feature_set, Scenario=$scenario"
            python "${SCRIPT_DIR}/test4.py" --probe "$probe" --feature_set "$feature_set" --scenario "$scenario" --output "results_all.csv"
        done
    done
done

echo "All runs completed."
