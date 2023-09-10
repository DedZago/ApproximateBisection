module load julia
# module load R
mkdir slurm-output
rm slurm-output/*
# julia -e "using Pkg; Pkg.add(\"DrWatson\"); using DrWatson; @quickactivate; Pkg.instantiate();"
FILE="cancel-submission.txt"
touch $FILE
INIT=1
END=100
for index in $(seq $INIT $END)
do
    export index
    echo seed is: $index
    job_id=$(sbatch --parsable run_job_parallel.slurm)
    if [[ $index -eq $INIT ]]
    then
        echo -n "scancel {$job_id.." >> $FILE 
    elif [[ $index -eq $END ]]
    then
        echo "$job_id}" >> $FILE 
    fi
done

