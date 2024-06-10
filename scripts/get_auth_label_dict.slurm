#!/bin/bash
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --mem 50G
#SBATCH -t 1-00:00
#SBATCH -o auth_label_dict.%N.%j.log
#SBATCH -e auth_label_dict.%N.%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_email_address
#SBATCH --account=def-*
##
#
python3 -u get_auth_label_dict.py -d /home/username/scratch/pdb_cif -o /home/username/scratch/interfacial_residues/files