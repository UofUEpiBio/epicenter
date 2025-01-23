#!/bin/sh
#SBATCH --job-name=ergms
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --time=00:05:00
#SBATCH --mem=5GB
#SBATCH --ntasks=2
#SBATCH --nodes=1

module load singularity
make singularity_render_chpc
