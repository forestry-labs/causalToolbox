#!/bin/bash
#SBATCH --cpus-per-task 1
R CMD BATCH --no-save Rand_TuneT.R Rand_TuneT_setup_1.Rout

