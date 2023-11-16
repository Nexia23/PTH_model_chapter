#!/bin/bash


#SBATCH --array=1-10 # how many tasks in the array
#SBATCH -c 1 # one CPU core per task
#SBATCH -o hellopatient59-%j-%a.out


srun python3 paramfitt_singlePatient.py $SLURM_ARRAY_TASK_ID 36


