module load r
module load gcc/4.8.5


for i in {1..11}; do sbatch BRC_Rand_TuneX.sh $i; done