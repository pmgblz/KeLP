#!/bin/bash
#
#SBATCH --job-name=ht1
#
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00

#save job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `

# run code
Rscript tune_gamma.R "whitenonbritish_unrelated" 'c("height")'  0.1 "TRUE" "FALSE"
Rscript tune_gamma.R "whitenonbritish_unrelated" 'c("height")'  0.2 "TRUE" "FALSE"
Rscript tune_gamma.R "whitenonbritish_unrelated" 'c("height")'  0.3 "TRUE" "FALSE"

#echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `



