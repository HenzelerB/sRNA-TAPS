#!/bin/bash
# =============================================================================
# sRNA_taps.bash — launch the sRNA-TAPS pipeline as a SLURM coordinator job.
#
# Run from your PROJECT directory (the one containing config.yaml):
#     bash /path/to/sRNA-TAPS/sRNA_taps.bash              # full pipeline
#     bash /path/to/sRNA-TAPS/sRNA_taps.bash --benchmark  # + rastair/asTair/Bismark
#
# This submits ONE small, long-walltime SLURM job that runs the Snakemake
# coordinator (via `srnataps run --slurm`). The coordinator then submits all
# compute steps as their own SLURM jobs using the bundled profile.
#
# Running the coordinator as a SLURM job (not on the login node) is deliberate:
# login-node processes get killed on disconnect / node maintenance, which
# silently aborts long runs. As a job, it survives until the pipeline finishes.
#
# Monitor:  squeue -u "$USER"
#           tail -f logs/coordinator_<jobid>.log
# =============================================================================
set -euo pipefail

# ── Parse optional --benchmark passthrough ───────────────────────────────────
BENCH_FLAG=""
for arg in "$@"; do
    case "$arg" in
        --benchmark) BENCH_FLAG="--benchmark" ;;
        *) echo "Unknown argument: $arg (only --benchmark is supported)"; exit 1 ;;
    esac
done

# ── Resolve paths ─────────────────────────────────────────────────────────────
PROJECT_DIR="$(pwd)"
CONFIG="${PROJECT_DIR}/config.yaml"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.yaml not found in current directory: ${PROJECT_DIR}"
    echo "Run this script from your project directory (the one holding config.yaml)."
    exit 1
fi
mkdir -p "${PROJECT_DIR}/logs"

echo "=============================================="
echo "sRNA-TAPS coordinator launcher"
echo "Project dir : ${PROJECT_DIR}"
echo "Config      : ${CONFIG}"
echo "Benchmark   : ${BENCH_FLAG:-(off)}"
echo "=============================================="

# ── Submit the coordinator as a SLURM job ─────────────────────────────────────
# Small footprint: the coordinator only orchestrates + waits on child jobs.
sbatch <<SBATCH_SCRIPT
#!/bin/bash
#SBATCH --job-name=srnataps_coord
#SBATCH --output=${PROJECT_DIR}/logs/coordinator_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=7-00:00:00
set -euo pipefail
cd "${PROJECT_DIR}"
srnataps run --slurm --configfile "${CONFIG}" ${BENCH_FLAG}
SBATCH_SCRIPT

echo "Coordinator job submitted. Track with: squeue -u \"\$USER\""
