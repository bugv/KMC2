#!/bin/bash

mkdir -p logs
mkdir -p results_first_processing/GA1_GB10_binary

for input_file in inputs/GA1_GB10_binary/*
do
    echo "Submitting job for $input_file"
    sbatch single_composition.slurm "$input_file"
done