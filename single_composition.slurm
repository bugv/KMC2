#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
#SBATCH --mem=400M
#SBATCH --output=logs/slurm-%A_%a.out

#SBATCH --array=0-5000%1000

input_file=$1  # passed when submitting
input_basename=$(basename "$input_file")  # strip path
input_name="${input_basename%.*}"          # strip extension

# Make output directory if not existing
mkdir -p results_first_processing/ideal_binary/${input_name}

# Run
python main.py -binary "$input_file" "results_first_processing/ideal_binary/${input_name}/${SLURM_ARRAY_TASK_ID}.h5"
