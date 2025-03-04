#!/bin/bash

################################################################################
# SLURM Job Parameters
################################################################################

#SBATCH --partition=sandbox            # Partition name
#SBATCH --time=0-04:00:00              # Format: DD-HH:MM:SS
#SBATCH --cpus-per-task=1              # Number of CPU cores per task
#SBATCH --mem=6400                     # Memory (in MB) per node

################################################################################
# Output and Error Files
################################################################################
# Use the "%j" token for job ID so each log file is uniquely named.

#SBATCH --error=logs/job_logs_%x.err
#SBATCH --output=logs/job_logs_%x.out

################################################################################
# Email Notifications (Uncomment & set up correct email if desired)
################################################################################

#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=kp779@hawaii.edu

################################################################################
# Environment Setup
################################################################################

# Declar session for logging
echo "This job name is $SLURM_JOB_NAME $SLURM_JOB_ID"

# Load your shell config (if needed) and activate conda environment
source ~/.zshrc
conda activate optoproject

################################################################################
# Pre-Processing: Pull latest GitHub
################################################################################


# Change directory to working directory 
cd ~/BRIDGE_Connection/

# Updating to current data
git pull
git checkout dev


################################################################################
# Main Execution
################################################################################

# Change directory to where your Python script is located
cd ~/BRIDGE_Connection/src/BRIDGE_Connection

# Run the script


python3 compliance.py cantilever ${SLURM_JOB_ID}




################################################################################
# Post-Processing: Commit and Push Logs to GitHub
################################################################################

# Change directory to working directory 
cd ~/BRIDGE_Connection/

# Make sure we’re on the logs branch (or create it if needed)
#git checkout logs || git checkout -b logs

# Copy any newly created log files into a subfolder or keep them in logs/
git add logs/job_logs_${SLURM_JOB_NAME}.*

# Commit changes (use a helpful commit message)
git commit -m "Add logs for job ${SLURM_JOB_NAME} ${SLURM_JOB_ID}"

# Push to remote (assumes SSH setup or a saved credential)
git push origin dev

echo "Logs pushed to 'logs' branch on GitHub..."
