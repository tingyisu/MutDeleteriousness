#!/bin/bash
#SBATCH -n 5
#SBATCH --mem-per-cpu 10G
#SBATCH -t 8:00:00
#SBATCH -o edgetic.%N.%j.log
#SBATCH -e edgetic.%N.%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
python -u get_edgetic_mutations.py -s /home/username/projects/def-*/username/modeller -c /home/username/scratch