'''
Writes commands to submit SLURM jobs to a BASH file
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse

# writes jobs to submit_residues_slurm_jobs.bash
def write_jobs(num_splits_for_hiunion, num_splits_for_intact):
	with open('submit_residues_slurm_jobs.bash', 'w') as f:
		f.write('#!/bin/bash' + '\n')
		# write hiunion SLURM jobs
		for i in range(1, int(num_splits_for_hiunion) + 1):
			f.write('sbatch ' + 'hiunion' + str(i) + '.slurm' + '\n')
		# write intact SLURM jobs
		for i in range(1, int(num_splits_for_intact) + 1):
			f.write('sbatch ' + 'intact' + str(i) + '.slurm' + '\n')


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-u', '--num_splits_for_hiunion') # number of slurm jobs that the HI-union data was split into
	parser.add_argument('-i', '--num_splits_for_intact') # number of slurm jobs that the IntAct data was split into
	args = parser.parse_args()

	write_jobs(args.num_splits_for_hiunion, args.num_splits_for_intact)

if __name__=='__main__':
	main()