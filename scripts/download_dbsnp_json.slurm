#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1 
#SBATCH --mem=50GB
#SBATCH -t 5-00:00
#SBATCH -o download_dbsnp_json.%N.%j.log
#SBATCH -e download_dbsnp_json.%N.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
##
#
# download dbsnp .json.bz files (one for each chromosome) first
chroms=('1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y')
for i in "${chroms[@]}"
do
	# if flatfile doesn't exist yet for chromosome i
	if [[ ! -f "refsnp-chr"$i".json.bz2" ]]
	then
		wget ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chr"$i".json.bz2
	fi
done