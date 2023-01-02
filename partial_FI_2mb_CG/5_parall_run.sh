#!/bin/sh
#
#
# Replace <ACCOUNT> with your account name before submitting.
#
#SBATCH --job-name=fi_cap    # The job name.
#SBATCH --gres=gpu:4
#SBATCH --array=1-5
#SBATCH -c 24
#SBATCH --time=1000:59:59              # The time the job will take to run.
#SBATCH --mem-per-cpu=16gb        # The memory the job will use per cpu core.
 
module load gcc-9.2.0
source /shared/apps/bin/GMXRC
python 5_parall_run.py $SLURM_ARRAY_TASK_ID

# End of script

