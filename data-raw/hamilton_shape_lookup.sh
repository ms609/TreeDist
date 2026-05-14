#!/bin/bash
#SBATCH --job-name=acid-shape
#SBATCH --partition=shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --array=10-12
#SBATCH --output=/nobackup/%u/acid-aphid/logs/%x_%A_%a.out
#SBATCH --error=/nobackup/%u/acid-aphid/logs/%x_%A_%a.err

set -euo pipefail
module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=/nobackup/$USER/acid-aphid/lib
export R_LIBS_USER=$LIB

n=$SLURM_ARRAY_TASK_ID

# nSim scales down as W(n) grows.  Targets: n=10 nSim=500 (~108 min), n=11
# nSim=200 (~3 h), n=12 nSim=100 (~7-8 h).  All comfortably under 24 h.
case $n in
  10) nSim=500 ;;
  11) nSim=200 ;;
  12) nSim=100 ;;
  *)  nSim=200 ;;
esac

out=/nobackup/$USER/acid-aphid/results/shape_lookup_n${n}.rds
Rscript /nobackup/$USER/acid-aphid/scripts/shape_lookup_n.R $n $nSim $out
