#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
#SBATCH --mem=400M
#SBATCH --output=logs/slurm-%A_%a.out

#SBATCH --array=0-499%100

input_file=$1  # passed when submitting
input_basename=$(basename "$input_file")  # strip path
input_name="${input_basename%.*}"          # strip extension

# Make output directory if not existing
mkdir -p results_first_processing/single_calculations/${input_name}

# Run 10 jobs sequentially
for i in $(seq 0 9); do
    task_id=$((SLURM_ARRAY_TASK_ID * 10 + i))
    python main.py -binary "$input_file" "results_first_processing/single_calculations/${input_name}/${task_id}.h5"
done
