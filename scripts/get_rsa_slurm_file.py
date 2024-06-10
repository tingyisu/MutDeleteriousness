'''
Processes and splits RSA (Relative Solvent Accessibility) jobs to SLURM files
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os
import os.path as osp
import argparse
import pickle
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump

class RSA:
	def __init__(self, script_dir, scratch_dir, account):
		self.account = account
		self.all_foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		self.script_dir = script_dir
		self.data_dir = osp.join(self.script_dir, 'files')
		# self.additional_foldx_buildmodel_mutations = pickle_load(osp.join(self.all_foldx_buildmodel_dir, 'additional_foldx_buildmodel_mutations.pickle'))
		self.split_folders = [f.path for f in os.scandir(self.all_foldx_buildmodel_dir) if f.is_dir() and 'split_' in f.path]
		self.additional_folders = [f.path for f in os.scandir(self.all_foldx_buildmodel_dir) if f.is_dir() and 'additional_' in f.path]
		# print(len(self.split_folders))
		# print(len(self.additional_folders))

	def split_jobs(self, folders, name, num_folders_per_job, pickle_file): # name = 'foldx_buildmodel_mutations' or 'foldx_buildmodel_mutations_additional'
		# get dictionary to keep track of mutations within each split folder
		split_dict = {} # key = split # (str), value = list of foldx_buildmodel_mutations in corresponding split
		# print(split_files, len(split_files))
		for file in folders:
			split_num = file.split('/')[-1].split('_')[-1].split('.')[0]
			split_dict[split_num] = pickle_load(osp.join(self.all_foldx_buildmodel_dir, name + '_' + split_num + '.pickle'))
		# print(len(split_files), len(self.foldx_buildmodel_mutations_split_dict['0']))
		pickle_dump(split_dict, pickle_file)
		num_parts = -1
		if len(folders) % num_folders_per_job == 0:
			num_parts = int(len(folders)/num_folders_per_job)
		else:
			num_parts = int(len(folders)/num_folders_per_job) + 1
		return split_dict, num_parts

	# writes RSA split/additional jobs to SLURM files
	def get_rsa_slurm_files(self, split_dict, name, num_folders_per_job, num_parts, num_days): # name = 'split' or 'additional'
		# combined 5 split jobs at a time
		split_num_list = list(split_dict.keys())
		# print(split_num_list)
		for i in range(num_parts):
			split_nums_for_part = []
			if i != num_parts-1:
				split_nums_for_part = split_num_list[i*num_folders_per_job:(i+1)*num_folders_per_job]
			else:
				split_nums_for_part = split_num_list[i*num_folders_per_job:len(split_num_list)]
			print(name, i, split_nums_for_part)

			with open(osp.join(self.script_dir, '_'.join(['calc_rsa', name, str(i) + '.slurm'])), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -N 1' + '\n')
				f.write('#SBATCH -n 5' + '\n')
				f.write('#SBATCH --mem 20G' + '\n')
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
				f.write('#SBATCH -o ' + '_'.join(['rsa', name, str(i) + '.%N.%j.log']) + '\n')
				f.write('#SBATCH -e ' + '_'.join(['rsa', name, str(i) + '.%N.%j.log']) + '\n')
				# get which split folder the mutation is in
				for split_num in split_nums_for_part:
					split_dir = osp.join(self.all_foldx_buildmodel_dir, name + '_' + split_num)
					for foldx_buildmodel_mutation in split_dict[split_num]:
						dirname = osp.join(split_dir, foldx_buildmodel_mutation)
						# get .pdb file with non-mutated residue, would be + '_1.pdb' for mutated residue (output of FoldX BuildModel)
						f.write('python calc_rsa.py ' + osp.join(dirname, foldx_buildmodel_mutation + '.pdb') + ' -o ' + osp.join(dirname, foldx_buildmodel_mutation)  + '\n')
						f.write('rm ' +  osp.join(dirname, foldx_buildmodel_mutation + '.asa.txt') + '\n')

	def get_rsa_slurm_files_all(self):
		num_folders_per_job = 5

		# for FoldX BuildModel mutations
		split_dict, num_parts = self.split_jobs(self.split_folders, 'foldx_buildmodel_mutations', num_folders_per_job, osp.join(self.all_foldx_buildmodel_dir, 'foldx_buildmodel_mutations_split_dict.pickle'))
		self.get_rsa_slurm_files(split_dict, 'split', num_folders_per_job, num_parts, '2')
		print('Number of RSA jobs for BuildModel:', num_parts)

		# for FoldX BuildModel additional mutations
		additional_dict, num_parts_additional = self.split_jobs(self.additional_folders, 'foldx_buildmodel_mutations_additional', num_folders_per_job, osp.join(self.all_foldx_buildmodel_dir, 'foldx_buildmodel_mutations_additional_dict.pickle'))
		self.get_rsa_slurm_files(additional_dict, 'additional', num_folders_per_job, num_parts_additional, '2')
		print('Number of RSA jobs for BuildModel additional:', num_parts_additional)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	parser.add_argument('-a', '--account') # def-*
	args = parser.parse_args()

	r = RSA(args.script_dir, args.scratch_dir, args.account)
	r.get_rsa_slurm_files_all()

if __name__=='__main__':
	main()