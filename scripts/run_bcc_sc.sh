#!/bin/bash

mkdir -p logs

# bcc
mkdir -p results_first_processing/gamma_spread_bcc #!!!

for input_file in inputs/gamma_spread_bcc/* #!!!
do
    input_basename=$(echo "$input_file" | awk -F'/' '{print $(NF-1)"/"$NF}')
    input_name="${input_basename%.*}"
    output_dir="results_first_processing/${input_name}"

    file_count=0
    if [ -d "$output_dir" ]; then
        file_count=$(find "$output_dir" -type f | wc -l)
    fi

    if [ "$file_count" -gt 5000 ]; then
        echo "Skipping $input_file (already has $file_count files)"
        continue
    fi

    echo "Submitting job for $input_file"
    sbatch scripts/single_composition.slurm "$input_file"
done

# sc
mkdir -p results_first_processing/gamma_spread_sc

for input_file in inputs/gamma_spread_sc/*
do
    input_basename=$(echo "$input_file" | awk -F'/' '{print $(NF-1)"/"$NF}')
    input_name="${input_basename%.*}"
    output_dir="results_first_processing/${input_name}"

    file_count=0
    if [ -d "$output_dir" ]; then
        file_count=$(find "$output_dir" -type f | wc -l)
    fi

    if [ "$file_count" -gt 5000 ]; then
        echo "Skipping $input_file (already has $file_count files)"
        continue
    fi

    echo "Submitting job for $input_file"
    sbatch scripts/single_composition.slurm "$input_file"
done

# unary 
mkdir -p results_first_processing/unary #!!!

for input_file in inputs/unary/* #!!!
do
    input_basename=$(echo "$input_file" | awk -F'/' '{print $(NF-1)"/"$NF}')
    input_name="${input_basename%.*}"
    output_dir="results_first_processing/${input_name}"

    file_count=0
    if [ -d "$output_dir" ]; then
        file_count=$(find "$output_dir" -type f | wc -l)
    fi

    if [ "$file_count" -gt 5000 ]; then
        echo "Skipping $input_file (already has $file_count files)"
        continue
    fi

    echo "Submitting job for $input_file"
    sbatch scripts/single_composition.slurm "$input_file"
done