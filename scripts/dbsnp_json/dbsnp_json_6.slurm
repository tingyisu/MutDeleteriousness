#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1 
#SBATCH --mem=50GB
#SBATCH -t 5-00:00
#SBATCH -o dbsnp_6.%N.%j.log
#SBATCH -e dbsnp_6.%N.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
##
#
# process dbsnp .json.bz files (one for each chromosome)
chroms=('21' '22' 'X' 'Y')
for i in "${chroms[@]}"
do
	python3 -u process_dbsnp_json.py -c $i -d /home/username/scratch/dbsnp/data
	# rm refsnp-chr"$i".json.bz2 # (optional to free up space)
done