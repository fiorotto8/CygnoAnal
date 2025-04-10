#!/bin/bash

# Arguments passed from Condor submission
param1=$1
param2=$2
param3=$3
param4=$4

# Run the first command
./INAF Reco/ArCF4_0deg.root . "$param1" "$param2" "$param3" "$param4"

# Run the second command
python fit_modulation.py --file Analysis_ArCF4_0deg.root