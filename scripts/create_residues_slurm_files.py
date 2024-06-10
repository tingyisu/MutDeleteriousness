'''
Creates SLURM files to find interfacial residues for PDB chain pairs
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse

# writes .slurm files to file
def create_slurm_jobs(interactome_name, curr_job_num, num_days, pdb_download_path, output_directory, email_address, account):
	name = interactome_name + str(curr_job_num)
	fname = name + '.slurm'
	python_script_name = interactome_name + '_residues_split.py'
	with open(fname, 'w') as f:
		f.write('#!/bin/bash' + '\n')
		f.write('#SBATCH -N 1' + '\n')
		f.write('#SBATCH -n 5' + '\n')
		f.write('#SBATCH --mem 10G' + '\n')
		f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n') # 10 days 
		f.write('#SBATCH -o ' + name + '.%N.%j.log' + '\n')
		f.write('#SBATCH -e ' + name + '.%N.%j.log' + '\n')
		f.write('#SBATCH --mail-type=END,FAIL' + '\n')
		f.write('#SBATCH --mail-user=' + str(email_address) + '\n')
		f.write('#SBATCH --account=' + account + '\n')
		f.write('##' + '\n')
		f.write('#' + '\n')
		f.write('python3 -u ' + python_script_name + ' -d ' + str(pdb_download_path) + ' -o ' + output_directory + ' -p ' + str(curr_job_num))

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--pdb_download_path') # e.g. /home/username/scratch/pdb_cif
	parser.add_argument('-o', '--output_directory') # e.g. /home/username/scratch/interfacial_residues/files
	parser.add_argument('-u', '--num_splits_for_hiunion') # number of slurm jobs to split HI-union data into
	parser.add_argument('-i', '--num_splits_for_intact') # number of slurm jobs to split IntAct data into
	parser.add_argument('-n', '--num_days') # number of days to run each slurm job
	parser.add_argument('-e', '--email_address') # your email address to send slurm job updates to
	parser.add_argument('-a', '--account') # account name/compute canada group (e.g. def-*)

	args = parser.parse_args()

	# HI-union SLURM jobs
	for i in range(1, int(args.num_splits_for_hiunion) + 1):
		create_slurm_jobs('hiunion', i, args.num_days, args.pdb_download_path, args.output_directory, args.email_address, args.account)

	# IntAct SLURM jobs
	for i in range(1, int(args.num_splits_for_intact) + 1):
		create_slurm_jobs('intact', i, args.num_days, args.pdb_download_path, args.output_directory, args.email_address, args.account)


if __name__=="__main__":
	main()