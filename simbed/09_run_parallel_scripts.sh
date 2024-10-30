#!/bin/bash

echo "Starting R scripts in parallel..."

# List of arguments
arguments=(1 2 3 4 5 6 7 8 9 10)

# Loop over arguments and run each R script in the background
for arg in "${arguments[@]}"
do
  Rscript 09_state_inference_parallel.R "$arg" &
done

# Wait for all background tasks to complete
wait

echo "All R scripts have finished running."
