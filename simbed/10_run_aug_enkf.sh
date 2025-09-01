#!/bin/bash

echo "Starting R scripts in parallel..."

# List of arguments
sample_inds=(1 3)

# Loop over arguments and run each R script in the background
for s in "${sample_inds[@]}"
do
    Rscript 10_enkf_stateaug.R "$s" &
done

# Wait for all background tasks to complete
wait

echo "All R scripts have finished running."
