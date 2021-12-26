#!/bin/bash
#SBATCH --job-name=array-job     # create a short name for your job
#SBATCH --output=slurm-%N.%j.out # STDOUT file
#SBATCH --error=slurm-%N.%j.err  # STDERR file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G         # memory per cpu-core (4G is default)
#SBATCH --time=03:50:00          # total run time limit (HH:MM:SS)
#SBATCH --array=0-700   # job array with index values 0, 1, 2, 3, 4
#SBATCH --mail-type=all          # send email on job start, end and fault
#SBATCH --mail-user=eivshina@princeton.edu


module load anaconda3/2020.11
source activate tt
python main.py $SLURM_ARRAY_TASK_ID
