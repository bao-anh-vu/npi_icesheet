#!/bin/bash

echo "Starting R scripts in parallel..."

# List of arguments
sample_inds=(7 8)

runs=(1 2 3 4 5)

# Loop over arguments and run each R script in the background
for s in "${sample_inds[@]}"
do
  for arg in "${runs[@]}"
  do
    Rscript 05_state_inference_parallel.R "$s" "$arg" &
  done
done

runs2=(6 7 8 9 10)

# Loop over arguments and run each R script in the background
for s in "${sample_inds[@]}"
do
  for arg in "${runs2[@]}"
  do
    Rscript 05_state_inference_parallel.R "$s" "$arg" &
  done
done


# Wait for all background tasks to complete
wait

echo "All R scripts have finished running."


# #!/bin/bash

# echo "Starting R scripts in parallel..."

# # List of arguments
# arguments=(1 2 3 4 5 6 7 8 9 10)

# # Loop over arguments and run each R script in the background
# for arg in "${arguments[@]}"
# do
#   Rscript 09_state_inference_parallel.R "$arg" &
# done

# # Wait for all background tasks to complete
# wait

# echo "All R scripts have finished running."
