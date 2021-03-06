#!/bin/bash

#SBATCH --job-name=6-obs
#SBATCH --error=/scratch/users/tballard/shum/giorgi.regions/reanalysis/6.err
#SBATCH --output=/scratch/users/tballard/shum/giorgi.regions/reanalysis/6.out
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=62GB
#SBATCH --mail-type=END
#SBATCH --mail-user=tballard@stanford.edu
#SBATCH --time=24:00:00
#SBATCH -p diffenbaugh

ml use /share/sw/modules/all
ml load R
ml load netCDF/4.3.0

Rscript 6-joint.frequency.1979.2014.R