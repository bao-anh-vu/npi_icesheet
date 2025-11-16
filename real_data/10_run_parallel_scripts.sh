#!/bin/bash

echo "Starting R scripts in parallel..."
sample_type="posterior"  # "prior" or "posterior"

# Create log file name dynamically based on sample_type
LOGFILE="run_log_${sample_type}.txt"
echo "Run started on $(date) for sample type: $sample_type" > "$LOGFILE"

# Parameters
start=1
end=50
group_size=10

FAILURES=()

run_group() {
    local samples=("$@")
    declare -A pids

    for s in "${samples[@]}"; do
        echo "Starting sample $s..."
        Rscript 10_state_inference_posterior.R "$sample_type" "$s" &
        pids[$!]=$s
    done

    for pid in "${!pids[@]}"; do
        wait "$pid"
        status=$?
        if [ $status -ne 0 ]; then
            echo "❌ Sample ${pids[$pid]} FAILED (exit code $status)" | tee -a "$LOGFILE"
            FAILURES+=("${pids[$pid]}")
        else
            echo "✅ Sample ${pids[$pid]} completed successfully" >> "$LOGFILE"
        fi
    done
}

# Build full list of samples
all_samples=($(seq $start $end))
num_samples=${#all_samples[@]}

group_index=1
for ((i=0; i<num_samples; i+=group_size)); do
    group=("${all_samples[@]:i:group_size}")
    echo "Running group $group_index: ${group[*]}"
    run_group "${group[@]}"
    echo "Group $group_index completed."
    ((group_index++))
done

echo "" >> "$LOGFILE"

if [ ${#FAILURES[@]} -eq 0 ]; then
    echo "🎉 All R scripts completed successfully!"
    echo "All runs successful." >> "$LOGFILE"
else
    echo "⚠ Some runs failed:"
    printf '%s\n' "${FAILURES[@]}"
    echo "Failed samples: ${FAILURES[*]}" >> "$LOGFILE"
fi

echo "Log available at $LOGFILE"
