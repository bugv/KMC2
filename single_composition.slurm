#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:03:00
#SBATCH --mem=400M
#SBATCH --output=logs/slurm-%A_%a.out

#SBATCH --array=0-5000%500

input_file=$1  # passed when submitting
input_basename=$(basename "$input_file")  # strip path
input_name="${input_basename%.*}"          # strip extension

# Make output directory if not existing

mkdir -p results_first_processing/GA1_GB10_binary/${input_name}

# Set output filename
output_filename="results_first_processing/ideal_binary/${input_name}/${SLURM_ARRAY_TASK_ID}.h5"

# Check if output file already exists
if [ -f "$output_filename" ]; then
    echo "Output file $output_filename already exists. Skipping."
    exit 0
fi

# Run
python main.py -binary "$input_file" "results_first_processing/GA1_GB10_binary/${input_name}/${SLURM_ARRAY_TASK_ID}.h5"
