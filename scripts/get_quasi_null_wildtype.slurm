#!/bin/bash
#SBATCH -n 5
#SBATCH --mem-per-cpu 10G
#SBATCH -t 2-00:00
#SBATCH -o quasi_null_wildtype.%N.%j.log
#SBATCH -e quasi_null_wildtype.%N.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
python -u get_quasi_null_wildtype.py -s /home/username/projects/def-*/username/modeller -c /home/username/scratch