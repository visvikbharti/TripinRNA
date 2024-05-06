
#!/bin/bash

#SBATCH --job-name=new_version_with_gap_minlen_7_job_05_12_23
#SBATCH --output=new_version_with_gap_minlen_7_job_05_12_23_%j.out
#SBATCH --cpus-per-task=18
#SBATCH --mem=32G
#SBATCH --time=32:00:00
#SBATCH --partition=compute

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR
echo "date=$(date)"
echo "Hostname=$(hostname -s)"
echo ""
echo "Number of Nodes Allocated=$SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated=$SLURM_NTASKS"
echo "Number of Cores/Tasks Allocated=$SLURM_CPUS_PER_TASK"
echo "Working Directory=$(pwd)"
echo "Working directory="$SLURM_SUBMIT_DIR


python src/with_gap.py
