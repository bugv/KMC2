#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:30:00
#SBATCH --mem=400M

#SBATCH --array=0-100

N=${SLURM_ARRAY_TASK_ID}

input_filename="inputs/composition_Al_0.50_Fe_0.50.json"
tmp_raw_kmc_output_file="tmp/results.json"
output_filename="results_first_processing/composition_Al_0.50_Fe_0.50/${SLURM_ARRAY_TASK_ID}.h5"

python main.py -binary $input_filename
python analysis_v2.py $tmp_raw_kmc_output_file $output_filename

