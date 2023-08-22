#!/bin/bash

#SBATCH --job-name=OPT-delayed
#SBATCH --mail-type=NONE                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dzago@ufl.edu     # Where to send mail	
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1GB                   # Job memory request
#SBATCH --account=qiu
#SBATCH --time=10-00:00:00               # Time limit hrs:min:sec
#SBATCH --qos=qiu
#SBATCH --output=slurm-output/%x_%A_%a.log   # Standard output and error log


# sleep 3m
# ./submit_generate_opt_data.sh
sleep 30m
./submit_optimization.sh