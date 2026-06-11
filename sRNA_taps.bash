#!/bin/bash
# =============================================================================
# sRNA_taps.bash
# Launch the sRNA-TAPS pipeline in a persistent tmux session.
#
# Usage (from the project directory, e.g. sRNA_TAPS_bio/):
#   bash /path/to/sRNA-TAPS/sRNA_taps.bash
#
# The Snakemake coordinator runs directly on the login node inside tmux
# and submits all compute jobs to SLURM via the profile/.
# This is stable — no SLURM timeout for the coordinator.
#
# Session management:
#   tmux ls                       — list sessions
#   tmux attach -t sRNA_TAPS      — reattach
#   tmux kill-session -t sRNA_TAPS — stop
# =============================================================================

SESSION_NAME="sRNA_TAPS"

# Activate conda environment
source /opt/miniforge3/etc/profile.d/conda.sh
source /opt/miniforge3/etc/profile.d/mamba.sh 2>/dev/null || true
conda activate RNA_taps

# Path to this script's directory (= sRNA-TAPS repo root)
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SNAKEFILE="${REPO_DIR}/workflow/Snakefile"
PROFILE="${REPO_DIR}/profile"

# Project directory is wherever the user runs this script from
PROJECT_DIR="$(pwd)"
CONFIG="${PROJECT_DIR}/config.yaml"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.yaml not found in current directory: ${PROJECT_DIR}"
    echo "Run this script from your project directory (sRNA_TAPS_bio/ or sRNA_TAPS_test/)"
    exit 1
fi

echo "=============================================="
echo "sRNA-TAPS pipeline launcher"
echo "Project dir : ${PROJECT_DIR}"
echo "Snakefile   : ${SNAKEFILE}"
echo "Profile     : ${PROFILE}"
echo "Config      : ${CONFIG}"
echo "tmux session: ${SESSION_NAME}"
echo "=============================================="

/usr/bin/tmux new-session -d -s "${SESSION_NAME}" \
    "source /opt/miniforge3/etc/profile.d/conda.sh && source /opt/miniforge3/etc/profile.d/mamba.sh 2>/dev/null; \
     conda activate RNA_taps && \
     cd ${PROJECT_DIR} && \
     /opt/apps/conda/bhenzeler/envs/RNA_taps/bin/snakemake \
         --snakefile ${SNAKEFILE} \
         --configfile ${CONFIG} \
         --profile ${PROFILE} \
         2>&1 | tee logs/snakemake_coordinator.log"

echo "Pipeline started in tmux session '${SESSION_NAME}'."
echo "Attach with: tmux attach -t ${SESSION_NAME}"
