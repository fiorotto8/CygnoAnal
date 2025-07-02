#!/bin/bash

# Path to directory with all outputs
BASEDIR="ParamsScan"
GLOBAL_OUTPUT="$BASEDIR/global_output.csv"
FILENAME="INAF_Ar_2ndreco.root" #name of original file

# Reset global output file and write CSV header
echo "run,modulation_factor_curve,error_curve,chi2_curve,modulation_factor_stokes,error_stokes,polarization_angle,polarization_angle_err,NPIP,w_facto,threshold,remove_noise_value" > "$GLOBAL_OUTPUT"

# Find all output directories, sort by run number
for output_dir in "$BASEDIR"/output_*; do
    # Extract run number (assuming output_DIR is output_XX)
    run_num=$(basename "$output_dir" | sed 's/output_//')
    echo "Analyzing $output_dir (run $run_num)"

    # Run the fit_modulation.py script, specify run number and output folder
    python3 fit_modulation.py --file "$output_dir/AnalysisMP_$FILENAME" --out_folder "$output_dir" --run "$run_num" --draw

    summary_file="$output_dir/output_${run_num}.txt"
    config_file="$BASEDIR/configFiles/config_${run_num}.txt"

    # Extract parameters from config file
    if [[ -f "$config_file" ]]; then
        NPIP=$(grep -E '^NPIP=' "$config_file" | cut -d'=' -f2)
        WFACTO=$(grep -E '^wfactor=' "$config_file" | cut -d'=' -f2)
        THRESHOLD=$(grep -E '^threshold=' "$config_file" | cut -d'=' -f2)
        REMOVENOISE=$(grep -E '^remove_noise_value=' "$config_file" | cut -d'=' -f2)
        echo "Parameters for run $run_num: NPIP=$NPIP, WFACTO=$WFACTO, THRESHOLD=$THRESHOLD, REMOVENOISE=$REMOVENOISE"
    else
        NPIP=""
        WFACTO=""
        THRESHOLD=""
        REMOVENOISE=""
        echo "WARNING: $config_file not found for run $run_num"
    fi

    # Check if summary file exists
    if [[ -f "$summary_file" ]]; then
        # Read the last line of the summary file
        read -r last_line < "$summary_file"
        # Extract values from the last line
        IFS=',' read -r run modulation_factor_curve error_curve chi2_curve modulation_factor_stokes error_stokes polarization_angle polarization_angle_err <<< "$last_line"
        # Append to global output CSV
        echo "$run_num,$modulation_factor_curve,$error_curve,$chi2_curve,$modulation_factor_stokes,$error_stokes,$polarization_angle,$polarization_angle_err,$NPIP,$WFACTO,$THRESHOLD,$REMOVENOISE" >> "$GLOBAL_OUTPUT"
    else
        echo "WARNING: $summary_file not found for run $run_num"
        # Append empty values for this run
        echo "$run_num,,," >> "$GLOBAL_OUTPUT"
    fi
    echo "Run $run_num analysis complete."
    echo "Results appended to $GLOBAL_OUTPUT"
    echo "----------------------------------------" 

done
