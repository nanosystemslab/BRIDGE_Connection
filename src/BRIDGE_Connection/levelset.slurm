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

#SBATCH --error=~/BRIDGE_Connection/logs/${JOB_NAME}_%j.err
#SBATCH --output=~/BRIDGE_Connection/logs/${JOB_NAME}_%j.out

################################################################################
# Email Notifications (Uncomment & set up correct email if desired)
################################################################################

#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
#SBATCH --mail-user=kp779@hawaii.edu

################################################################################
# Environment Setup
################################################################################

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
# python3 compliance.py cantilever ${SLURM_JOB_ID}
python3 compliance.py half_wheel ${SLURM_JOB_ID}


################################################################################
# Post-Processing: Commit and Push Logs to GitHub
################################################################################

# Change directory to working directory 
cd ~/BRIDGE_Connection/

# Make sure we’re on the logs branch (or create it if needed)
git checkout logs || git checkout -b logs

# Copy any newly created log files into a subfolder or keep them in logs/
# For example, if your logs are already in `~/BRIDGE_Connection/logs`,
# you can just commit them directly:
git add logs/${JOB_NAME}_${SLURM_JOB_ID}.*

# Commit changes (use a helpful commit message)
git commit -m "Add logs for job ${SLURM_JOB_ID} (${JOB_NAME})"

# Push to remote (assumes SSH setup or a saved credential)
git push origin logs

echo "Logs pushed to 'logs' branch on GitHub..."
