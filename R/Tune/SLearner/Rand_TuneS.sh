#!/bin/bash
#SBATCH --cpus-per-task 24
R CMD BATCH --no-save Rand_TuneS.R Rand_TuneS_setup_9.Rout

