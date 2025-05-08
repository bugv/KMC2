#!/bin/bash

mkdir -p logs
mkdir -p results_first_processing/gamma_spread

for input_file in inputs/gamma_spread/*
do
    echo "Submitting job for $input_file"
    sbatch scripts/single_composition.slurm "$input_file"
done