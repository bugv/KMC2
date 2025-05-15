#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
#SBATCH --mem=400M
#SBATCH --output=logs/slurm-%A_%a.out

#SBATCH --array=0-499%100

input_file=$1  # passed when submitting
input_basename=$(echo "$input_file" | awk -F'/' '{print $(NF-1)"/"$NF}')  # extract last two path components
input_name="${input_basename%.*}"          # strip extension

# Make output directory if not existing

mkdir -p results_first_processing/${input_name}

# Set output filename
output_filename="results_first_processing/${input_name}/${SLURM_ARRAY_TASK_ID}.h5"

# Check if output file already exists
if [ -f "$output_filename" ]; then
    echo "Output file $output_filename already exists. Skipping."
    exit 0
fi

# Run 10 jobs sequentially
for i in $(seq 0 9); do
    task_id=$((SLURM_ARRAY_TASK_ID * 10 + i))
    python main.py -binary "$input_file" "results_first_processing/${input_name}/${task_id}.h5"
done
