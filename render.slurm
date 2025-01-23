#!/bin/sh
#SBATCH --job-name=ergms
#SBATCH --output=ergms.out
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --time=48:00:00
#SBATCH --mem=10GB
#SBATCH --ntasks=4
#SBATCH --nodes=1
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --mail-type=ALL

module load singularity
make singularity_render_chpc
