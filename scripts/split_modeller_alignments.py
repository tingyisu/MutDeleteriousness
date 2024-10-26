'''
Split MODELLER alignments and pickle them to separate files/jobs
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os
import os.path as osp
import argparse
import pickle
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump

class Alignment:
	def __init__(self, script_dir, pdb_download_dir, email, account):
		self.script_dir = script_dir
		self.data_dir = osp.join(script_dir, 'files')
		self.alignments_dir = osp.join(script_dir, 'alignments')
		check_create_dir(self.alignments_dir)
		self.results_dir = osp.join(script_dir, 'results')
		check_create_dir(self.results_dir)
		self.pdb_download_dir = pdb_download_dir
		self.email = email
		self.account = account
		self.blast_alignments = pickle_load(osp.join(self.data_dir, 'converted_blast_alignments.pickle')) 
		self.num_alignments_per_job = 1000
		self.num_jobs = 0

	def create_uniprot_ali_file(self): # put uniprot sequence into .ali format specified by MODELLER
		for (uniprot, pdb_chain) in self.blast_alignments:
			pdb, chain = pdb_chain.split('_')
			fname = '_'.join([uniprot, pdb_chain])
			qstart, qend, auth_sstart, auth_send, qalignment, converted_salignment = self.blast_alignments[(uniprot, pdb_chain)]
			print('Creating .ali file for', uniprot, pdb_chain + '...')
			fname = osp.join(self.alignments_dir, fname + '.ali')
			cif_fname = osp.join(self.pdb_download_dir, pdb + '.cif')
			# if not check_exist(fname):
			with open(fname, 'w') as f:
				f.write('>P1;' + pdb_chain + '\n')
				f.write('structure:' + cif_fname + ':' + auth_sstart + ':' + chain + ':' + auth_send + ':' + chain + '::::' + '\n')
				f.write(converted_salignment + '*' + '\n' + '\n')
				f.write('>P1;' + uniprot + '\n')
				f.write('sequence:' + uniprot + ':' + qstart + '::' + qend + ':::::' '\n')
				f.write(qalignment + '*' + '\n')

	def split_blast_alignments(self): # splits self.blast_alignments into different parts (for running on compute canada) dictionary of blast alignments
		pairs_not_looked_at_yet = []
		for (uniprot, pdb_chain) in self.blast_alignments:
			fname = '_'.join([uniprot, pdb_chain])
			if not check_exist(osp.join(self.results_dir, fname + '.B99990001.pdb')):
				pairs_not_looked_at_yet.append((uniprot, pdb_chain))
		num_parts = int(len(pairs_not_looked_at_yet)/self.num_alignments_per_job) + 1
		print('Number of alignments to model:', len(pairs_not_looked_at_yet))
		self.num_jobs = num_parts
		total_num = 0
		total = []
		for i in range(num_parts):
			alignments_for_part = []
			if i != num_parts-1:
				alignments_for_part = pairs_not_looked_at_yet[i*self.num_alignments_per_job:(i+1)*self.num_alignments_per_job]
			else:
				alignments_for_part = pairs_not_looked_at_yet[i*self.num_alignments_per_job:len(pairs_not_looked_at_yet)]

			# pickle
			pickle_dump(alignments_for_part, osp.join(self.alignments_dir, 'blast_alignments_for_modeller_' + str(i) + '.pickle'))

			print('Part', i, 'contains', len(alignments_for_part), 'BLAST alignment pairs')
			total_num += len(alignments_for_part)
			total += alignments_for_part

		print('Total number of pairs:', total_num)
		if sorted(total) != sorted(pairs_not_looked_at_yet):
			print('ERROR! Something went wrong when splitting the BLAST alignments up!')

	def create_slurm_jobs(self):
		for i in range(self.num_jobs):
			with open(osp.join(self.script_dir, 'run_modeller_part_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -n 5' + '\n')
				f.write('#SBATCH --mem-per-cpu 10G' + '\n')
				f.write('#SBATCH -t 15:00:00' + '\n') # 15 hours
				f.write('#SBATCH -o run_modeller_part_' + str(i) + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e run_modeller_part_' + str(i) + '.%N.%j.log' + '\n')
				f.write('#SBATCH --mail-type=END,FAIL' + '\n') 
				f.write('#SBATCH --mail-user=' + self.email + '\n')
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('DIR=' + self.script_dir + '\n')
				f.write('python -u run_modeller_all.py -s "$DIR" -n ' + str(i) + '\n')


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  # /home/username/projects/def-*/username/modeller
	parser.add_argument('-p', '--pdb_download_dir') # /home/username/scratch/pdb_cif
	parser.add_argument('-e', '--email') # your email address
	parser.add_argument('-a', '--account') # def-* (compute canada project allocation/group name)
	args = parser.parse_args()

	a = Alignment(args.script_dir, args.pdb_download_dir, args.email, args.account)
	a.create_uniprot_ali_file()
	a.split_blast_alignments()
	a.create_slurm_jobs()


if __name__=='__main__':
	main()
