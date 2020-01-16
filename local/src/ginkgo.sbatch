#!/bin/bash

#SBATCH --job-name=ginkgo
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marilisa.montemurro@polito.it
#SBATCH --partition=global
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=ginkgo_%j.log
#SBATCH --mem-per-cpu=4096M


python /home/mmontemurro/PhyliCS/local/src/phylics.py --run_cnv  --beds_dir ../../dataset/10x_breast/breast_tissue_A-sc/ --genome GrCh38 --binning variable_500000_101_bwa --init_ginkgo
