#!/bin/bash

#SBATCH --job-name=RSADA_PARALLEL
#SBATCH --mail-type=NONE                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=dzago@ufl.edu     # Where to send mail	
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=5 
#SBATCH --hint=multithread             # Allow the use of in-core multi-threading
#SBATCH --mem=10GB                   # Job memory request
#SBATCH --time=95:00:00               # Time limit hrs:min:sec
#SBATCH --account=qiu
#SBATCH --qos=qiu-b
#SBATCH --output=slurm-output/%x_%A_%a.log   # Standard output and error log

julia -t 5 scripts/run-parallel.jl --index $index

