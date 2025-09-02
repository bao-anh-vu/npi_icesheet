#!/bin/bash

echo "Starting R scripts in parallel..."

# List of arguments
sample_inds=(1 6 7 8 9 10)

# Loop over arguments and run each R script in the background
for s in "${sample_inds[@]}"
do
    # Redirect stdout and stderr to a log file per sample
    Rscript 10_enkf_stateaug.R "$s" > "log/aug_enkf_log_s${s}.log" 2>&1 &
done

# Wait for all background tasks to complete
wait

echo "All R scripts have finished running."
