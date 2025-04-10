ssh 

    ssh bugnion@helvetios.hpc.epfl.ch

    cd /scratch/bugnion

DO NOT FORGET TO ACTIVATE VENV

    source venvs/kmc2/bin/activate


slurm things 

    Squeue 

    scancel -u $USER

    scancel -u $USER -t PENDING

    
number of files in each subdir of cwd 

    for d in */ ; do   echo "$d: $(find "$d" -type f | wc -l) files"; done