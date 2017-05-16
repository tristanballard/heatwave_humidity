#!/bin/bash

# This script runs multiple R scripts in parallel. Each R script (e.g. tmax.lon.loop.1.R)
# is applying a function (quantile regression) to every pixel in a 5268x4823 tmax grid
# from the Daymet tmax values that are daily from 1980-2014. This script runs the regression
# one month at a time. The first R script tells R to apply the function to, for a fixed 
# column (longitude), every row (latitude). Each R script is set to loop over 200 lon
# values each, leading to 24 R scripts. This scripts runs 5 of those 24 in parallel. 

#SBATCH --job-name=regrid
#SBATCH --error=/scratch/users/tballard/errors.1.err
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=100GB
#SBATCH --mail-type=END
#SBATCH --mail-user=tballard@stanford.edu
#SBATCH --time=12:00:00
#SBATCH -p owners

ml use /share/sw/modules/all
ml load R
ml load udunits
ml load netCDF/4.3.0

Rscript convert.to.rds.1.R

# parallel -j0 Rscript :::: <(ls tmax.lon.loop.{1..5}.R)
                           
# parallel tells the machine to run in parallel..
# -j0 tells it to run as many jobs as possible
# Rscript is the terminal command for running a .R file