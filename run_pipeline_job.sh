#!/bin/bash
#SBATCH --job-name=rnaseq_pipeline
#SBATCH --partition=extended-40core
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# Script directory (where sbatch was submitted from)
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"

# Create logs directory if it doesn't exist
mkdir -p "${SCRIPT_DIR}/logs"

# Print job info
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURMD_NODENAME"
echo "Start Time: $(date)"
echo "Script Directory: ${SCRIPT_DIR}"
echo "=========================================="

# Load Conda
source ~/miniconda3/etc/profile.d/conda.sh

# Activate R environment
conda activate rnaseq-r-env

# Change to script directory
cd "${SCRIPT_DIR}"

# Run the pipeline (choose one):
# Option 1: Run with padj only (default)
# Rscript run_pipeline.R 2>&1 | tee logs/pipeline_run_${SLURM_JOB_ID}.log

# Option 2: Run with p-value only
# Rscript run_pipeline.R --use-pvalue 2>&1 | tee logs/pipeline_run_${SLURM_JOB_ID}.log

# Option 3: Run both padj and p-value
Rscript run_pipeline.R --both 2>&1 | tee logs/pipeline_run_${SLURM_JOB_ID}.log

# Print completion info
echo "=========================================="
echo "End Time: $(date)"
echo "Exit Code: $?"
echo "=========================================="
