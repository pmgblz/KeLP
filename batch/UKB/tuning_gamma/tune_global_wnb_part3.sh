#!/bin/bash
#
#SBATCH --job-name=gt3
#
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00

#save job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `

# run code
Rscript tune_gamma.R "whitenonbritish_unrelated" 'c("platelet", "platelet_volume", "platelet_width", "platelet_crit")'  0.7 "TRUE" "FALSE"
Rscript tune_gamma.R "whitenonbritish_unrelated" 'c("platelet", "platelet_volume", "platelet_width", "platelet_crit")'  0.8 "TRUE" "FALSE"


#echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `



