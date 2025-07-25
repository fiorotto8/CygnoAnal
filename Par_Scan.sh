#!/bin/bash

# Define start, stop, and step for each parameter
NPIP_start=0; NPIP_stop=100; NPIP_step=10
wfactor_start=0; wfactor_stop=6; wfactor_step=1
threshold_start=0; threshold_stop=5; threshold_step=1
remove_noise_start=0; remove_noise_stop=10; remove_noise_step=1

# Calculate the number of combinations
num_NPIP=$(( (NPIP_stop - NPIP_start) / NPIP_step + 1 ))
num_wfactor=$(( (wfactor_stop - wfactor_start) / wfactor_step + 1 ))
num_threshold=$(( (threshold_stop - threshold_start) / threshold_step + 1 ))
num_noise=$(( (remove_noise_stop - remove_noise_start) / remove_noise_step + 1 ))
total_combinations=$(( num_NPIP * num_wfactor * num_threshold * num_noise ))

echo "Total combinations to be generated: $total_combinations"
read -p "Press Enter to continue..."

# Generate ranges automatically and shuffle the values
NPIP_values=($(seq $NPIP_start $NPIP_step $NPIP_stop | shuf))
wfactor_values=($(seq $wfactor_start $wfactor_step $wfactor_stop | shuf))
threshold_values=($(seq $threshold_start $threshold_step $threshold_stop | shuf))
remove_noise_values=($(seq $remove_noise_start $remove_noise_step $remove_noise_stop | shuf))

# Path to the base (template) config file
BASE_CONFIG="directionalityMPConfig_template.txt"
BASEDIR="ParamsScan"

rm -rf "$BASEDIR"
# Ensure config and output base directories exist
mkdir -p "$BASEDIR/configFiles"

runNum=0
run_skip=0  # Skip the first 6500 runs
# Loop over all combinations
for NPIP in "${NPIP_values[@]}"; do
    for wfactor in "${wfactor_values[@]}"; do
        for threshold in "${threshold_values[@]}"; do
            for noise in "${remove_noise_values[@]}"; do

                if (( runNum < run_skip )); then
                    runNum=$((runNum + 1))
                    continue
                fi

                # 1. Create a new config file for this parameter set
                config_file="$BASEDIR/configFiles/config_${runNum}.txt"
                cp "$BASE_CONFIG" "$config_file"
                sed -i "s/^NPIP=.*/NPIP=${NPIP}/"               "$config_file"
                sed -i "s/^wfactor=.*/wfactor=${wfactor}/"       "$config_file"
                sed -i "s/^threshold=.*/threshold=${threshold}/" "$config_file"
                sed -i "s/^remove_noise_value=.*/remove_noise_value=${noise}/" "$config_file"

                # 2. Set a unique output directory for this run and ensure it exists
                output_dir="$BASEDIR/output_${runNum}"
                mkdir -p "$output_dir"
                sed -i "s|^output_dir=.*|output_dir=${output_dir}|" "$config_file"

                # 3. Run the analysis program with the new config
                echo "Running NPIP=${NPIP}, wfactor=${wfactor}, threshold=${threshold}, noise=${noise} (runNum=${runNum})..."
                ./INAF_MP "$config_file"

                runNum=$((runNum + 1))
            done
        done
    done
done
done