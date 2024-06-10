#!/bin/bash
#SBATCH -n 5
#SBATCH --mem-per-cpu 10G
#SBATCH -t 3-00:00
#SBATCH -o prepare_modeller_pdb.%N.%j.log
#SBATCH -e prepare_modeller_pdb.%N.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
DIR=/home/username/projects/def-*/username/modeller
python -u prepare_modeller_pdb_files_for_foldx.py -s $DIR -c /home/username/scratch -a def-* -f foldx_20231231
