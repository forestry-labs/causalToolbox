#!/bin/bash
# Job name:
#SBATCH --job-name=TuneSX
#
# Account:
#SBATCH --account=fc_praxis
#
# Partition:
#SBATCH --partition=savio2
#
# Request one node:
#SBATCH --nodes=1
#
# Specify number of tasks for use case (example):
#SBATCH --ntasks-per-node=1
#
# Processors per task:
#SBATCH --cpus-per-task=24
#
# Wall clock limit:
#SBATCH --time=72:00:00
#
## Command(s) to run (example):#!/bin/bash
R CMD BATCH --no-save --args -$1 Rand_TuneX.R Rand_TuneX_setup_$1.Rout


