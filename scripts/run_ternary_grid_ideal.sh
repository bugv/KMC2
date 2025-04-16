#!/bin/bash

mkdir -p logs
mkdir -p results_first_processing/ideal_ternary

for input_file in inputs/ideal_ternary/*
do
    echo "Submitting job for $input_file"
    sbatch scripts/single_composition.slurm "$input_file"
done