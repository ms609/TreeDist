#!/bin/bash
#SBATCH --job-name=acid-build
#SBATCH --partition=shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=/nobackup/%u/acid-aphid/logs/%x_%j.out
#SBATCH --error=/nobackup/%u/acid-aphid/logs/%x_%j.err

set -euo pipefail
module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=/nobackup/$USER/acid-aphid/lib
export R_LIBS_USER=$LIB

cd /nobackup/$USER/acid-aphid/repos/TreeTools
git pull origin acid-aphid
rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library=$LIB TreeTools_*.tar.gz
rm -f TreeTools_*.tar.gz

cd /nobackup/$USER/acid-aphid/repos/TreeDist
git pull origin acid-aphid
rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library=$LIB TreeDist_*.tar.gz
rm -f TreeDist_*.tar.gz

echo "Build complete."
Rscript -e ".libPaths(c('$LIB', .libPaths())); cat('TreeTools', as.character(packageVersion('TreeTools')), '\n'); cat('TreeDist', as.character(packageVersion('TreeDist')), '\n')"
