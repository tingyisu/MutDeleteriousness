#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1 
#SBATCH --mem=50GB
#SBATCH -t 3-00:00
#SBATCH -o dbsnp_5.%N.%j.log
#SBATCH -e dbsnp_5.%N.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
##
#
# process dbsnp .json.bz files (one for each chromosome)
chroms=('17' '18' '19' '20')
for i in "${chroms[@]}"
do
	python3 -u process_dbsnp_json.py -c $i -d /home/username/scratch/dbsnp/data
	# rm refsnp-chr"$i".json.bz2 # (optional to free up space)
done