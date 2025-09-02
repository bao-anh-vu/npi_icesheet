#!/bin/bash

echo "Comparing samples..."

# List of arguments
arguments=(7 8)

# Loop over arguments and run each R script in the background
for arg in "${arguments[@]}"
do
  Rscript 11_compare_results.R "$arg" &
done

# Wait for all background tasks to complete
wait

echo "All R scripts have finished running."
