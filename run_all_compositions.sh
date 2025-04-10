#!/bin/bash

mkdir -p logs

for input_file in inputs/*
do
    echo "Submitting job for $input_file"
    sbatch single_composition.slurm "$input_file"
done