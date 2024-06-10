#!/bin/bash
#SBATCH -n 5
#SBATCH --mem-per-cpu 10G
#SBATCH -t 1-00:00
#SBATCH -o additional_foldx.%N.%j.log
#SBATCH -e additional_foldx.%N.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
DIR=/home/username/projects/def-*/username/modeller
python -u get_additional_foldx_buildmodel_mutations.py -s $DIR -c /home/username/scratch -a def-* -f foldx_20231231
