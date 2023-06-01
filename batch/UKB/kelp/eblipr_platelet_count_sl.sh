#!/bin/bash
#
#SBATCH --job-name=ep
#
#SBATCH --mem-per-cpu=120G
#SBATCH --partition=owners,normal,candes
#SBATCH --time=06:00:00

#save job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `

# run code
Rscript eblipr_UKB.R "british_unrelated" 'c("platelet")' 0.1 "alpha / 8" "FALSE" "FALSE" "FALSE" 


#echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `


