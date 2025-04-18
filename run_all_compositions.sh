#!/bin/bash

mkdir -p logs
mkdir -p results_first_processing/GA1_GB100_binary

for input_file in inputs/GA1_GB100_binary/*
do
    echo "Submitting job for $input_file"
    sbatch --nice=10000 single_composition.slurm "$input_file"
done