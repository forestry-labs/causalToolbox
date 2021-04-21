#!/bin/bash
#SBATCH --cpus-per-task 2
R CMD BATCH --no-save --args -$1 Rand_TuneX.R Rand_TuneX_setup_$1.Rout

