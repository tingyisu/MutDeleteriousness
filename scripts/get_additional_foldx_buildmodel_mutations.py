'''
Processing additional mutations (IR mutations that are non-edgetic) to run FoldX BuildModel on
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os
import os.path as osp
import argparse
import pickle
from Bio.PDB import PDBIO, PDBParser
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump, get_missense_info, write_missense_info_to_file

class AdditionalMutations:
	def __init__(self, script_dir, scratch_dir, account, foldx_executable_name):
		self.account = account
		self.foldx_executable_name = foldx_executable_name
		self.foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		self.script_dir = script_dir
		self.additional_foldx_buildmodel_mutations = pickle_load(osp.join(self.foldx_buildmodel_dir, 'additional_foldx_buildmodel_mutations.pickle'))
		self.uniprot_alignment_additional_dict = pickle_load(osp.join(self.script_dir, 'files', 'uniprot_alignment_dict_additional.pickle')) # key = (uniprot, pdb_chain), value = dict of modeller residue numbering (starting from one), and original uniprot residue numbering key-value pairs
		self.uniprot_alignment_dict = pickle_load(osp.join(self.script_dir, 'files', 'uniprot_alignment_dict.pickle'))
		self.uniprot_alignment_dict.update(self.uniprot_alignment_additional_dict) # merge
		self.num_mutations_per_job = 2000
		self.foldx_buildmodel_additional_num_splits = int(len(self.additional_foldx_buildmodel_mutations)/2000) + 1
		# self.foldx_buildmodel_additional_dir = osp.join(foldx_buildmodel_dir, 'additional')
		# check_create_dir(self.foldx_buildmodel_additional_dir)
		self.results_dir = osp.join(script_dir, 'results')

	# need to do this in order to avoid value error with overlapping residue numbers (https://github.com/biopython/biopython/issues/1551)
	# creates a residue_id_dict to map modeller numbered residues to original Uniprot numbered residues
	def get_residue_id_dict(self, model, residue_id_list):
		residue_id_dict = {}
		for residue_id in residue_id_list:
			list_residue_id = list(residue_id)
			original_res_num = list_residue_id[1]
			actual_res_num = self.uniprot_alignment_dict[model][original_res_num]
			# .pdb formats cannot handle cases where mut_pos_in_uniprot_protein > 9999, so do the following
			if actual_res_num > 10000:
				list_residue_id[1] = actual_res_num % 10000
			else:
				list_residue_id[1] = actual_res_num
			new_residue_id = tuple(list_residue_id)
			residue_id_dict[residue_id] = new_residue_id
		return residue_id_dict

	def create_parts(self, num_splits, dirname, mutations, name):
		total_num = 0
		total = []
		for i in range(num_splits):
			check_create_dir(osp.join(dirname, 'additional_' + str(i)))
			mutations_for_part = []
			if i != num_splits:
				mutations_for_part = mutations[i*self.num_mutations_per_job:(i+1)*self.num_mutations_per_job]
			else:
				mutations_for_part = mutations[i*self.num_mutations_per_job:len(mutations)]

			# pickle
			pickle_dump(mutations_for_part, osp.join(dirname, name + '_' + str(i) + '.pickle'))

			print('Part', i, 'contains', len(mutations_for_part), name)
			total_num += len(mutations_for_part)
			total += mutations_for_part

		print('Total number of mutations:', total_num)
		if sorted(total) != sorted(mutations):
			print('ERROR! Something went wrong when splitting the mutations!')

	def split_foldx_mutations(self):
		# first check whether or not a split exists already
		if not check_exist(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_additional_0.pickle')):
			self.create_parts(self.foldx_buildmodel_additional_num_splits , self.foldx_buildmodel_dir, self.additional_foldx_buildmodel_mutations, 'foldx_buildmodel_mutations_additional')

	def create_foldx_buildmodel_pdb_files(self):
		parser = PDBParser(QUIET=True)
		io = PDBIO()
		for i in range(self.foldx_buildmodel_additional_num_splits):
			print('Creating PDB files for FoldX BuildModel additional mutations part:', i)
			foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_additional_' + str(i) + '.pickle'))
			num_mutations = 1
			for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
				model, mutation = foldx_buildmodel_mutation.split('-')
				model_pdb_file = osp.join(self.results_dir, model + '.B99990001.pdb')
				dirname = osp.join(self.foldx_buildmodel_dir, 'additional_' + str(i), foldx_buildmodel_mutation)
				# create new directory for each foldx_buildmodel_mutation
				check_create_dir(dirname)
				new_pdb_file = osp.join(dirname, foldx_buildmodel_mutation + '.pdb')
				
				# get model1 and model2 structures
				structure = parser.get_structure(model, model_pdb_file)

				# rename residue numbers in model1 to uniprot numbering
				# get residue list for model1
				residue_id_list = []
				for residue in structure[0]['A']:
					residue_id_list.append(residue.id)

				# get residue id dict for model1
				# need to do this in order to not get value error with overlapping residue numbers (https://github.com/biopython/biopython/issues/1551)
				residue_id_dict = self.get_residue_id_dict(model, residue_id_list)
				for residue_id in residue_id_list[::-1]: # change residue number in reverse to avoid overlapping residue numbers
					structure[0]['A'][residue_id].id = residue_id_dict[residue_id]

				io.set_structure(structure)
				io.save(new_pdb_file)
					
				print('Done:', str(num_mutations) + '/' + str(len(foldx_buildmodel_mutations_for_part)))
				num_mutations += 1

	def create_foldx_buildmodel_config_files(self):
		for i in range(self.foldx_buildmodel_additional_num_splits):
			print('Creating FoldX BuildModel additional configuration files for part:', i)
			foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_additional_' + str(i) + '.pickle'))
			for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
				dirname = osp.join(self.foldx_buildmodel_dir, 'additional_' + str(i), foldx_buildmodel_mutation)
				_, foldx_mutation = foldx_buildmodel_mutation.split('-')
				# create mutant file
				with open(osp.join(dirname, 'individual_list.txt'), 'w') as f_mut:
					f_mut.write(foldx_mutation + ';')
				# write BuildModel config file
				with open(osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg'), 'w') as f:
					f.write('command=BuildModel' + '\n')
					f.write('pdb-dir=' + dirname + '\n')
					f.write('output-dir=' + dirname + '\n')
					f.write('pdb=' + foldx_buildmodel_mutation + '.pdb' + '\n')
					f.write('mutant-file=' + osp.join(dirname, 'individual_list.txt') + '\n')
					# the rest of the parameters are all default

	def create_foldx_buildmodel_slurm_job(self, num_days):
		for i in range(self.foldx_buildmodel_additional_num_splits):
			print('Creating FoldX BuildModel SLURM job for part:', i)
			with open(osp.join(self.foldx_buildmodel_dir, 'additional_' + str(i), 'foldx_buildmodel_additional_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -N 1' + '\n')
				f.write('#SBATCH -n 5' + '\n')
				f.write('#SBATCH --mem 20G' + '\n')
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
				f.write('#SBATCH -o ' + 'foldx_buildmodel_additional' + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e ' + 'foldx_buildmodel_additional' + '.%N.%j.log' + '\n')
				foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_additional_' + str(i) + '.pickle'))
				for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
					dirname = osp.join(self.foldx_buildmodel_dir, 'additional_' + str(i), foldx_buildmodel_mutation)
					f.write('./' + self.foldx_executable_name + ' -f ' + osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg') + '\n')
					# remove all intermediate/uneeded output files to conserve space (# of files) in scratch directory (maximum # of files = 1000k)
					f.write('rm ' + osp.join(dirname, 'Average_' + foldx_buildmodel_mutation + '.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'PdbList_' + foldx_buildmodel_mutation + '.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, foldx_buildmodel_mutation + '_1.pdb') + '\n')
					f.write('rm ' + osp.join(dirname, 'Raw_' + foldx_buildmodel_mutation + '.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'WT_' + foldx_buildmodel_mutation + '_1.pdb') + '\n')

	def create_foldx_buildmodel_config_files(self):
		for foldx_buildmodel_mutation in self.additional_foldx_buildmodel_mutations:
			dirname = osp.join(self.foldx_buildmodel_additional_dir, foldx_buildmodel_mutation)
			_, foldx_mutation = foldx_buildmodel_mutation.split('-')
			# create mutant file
			with open(osp.join(dirname, 'individual_list.txt'), 'w') as f_mut:
				f_mut.write(foldx_mutation + ';')
			# write BuildModel config file
			with open(osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg'), 'w') as f:
				f.write('command=BuildModel' + '\n')
				f.write('pdb-dir=' + dirname + '\n')
				f.write('output-dir=' + dirname + '\n')
				f.write('pdb=' + foldx_buildmodel_mutation + '.pdb' + '\n')
				f.write('mutant-file=' + osp.join(dirname, 'individual_list.txt') + '\n')
				# the rest of the parameters are all default

	def create_foldx_buildmodel_slurm_job(self, num_days):
		with open(osp.join(self.foldx_buildmodel_additional_dir, 'foldx_buildmodel_additional.slurm'), 'w') as f:
			f.write('#!/bin/bash' + '\n')
			f.write('#SBATCH -N 1' + '\n')
			f.write('#SBATCH -n 5' + '\n')
			f.write('#SBATCH --mem 20G' + '\n')
			f.write('#SBATCH --account=' + self.account + '\n')
			f.write('#SBATCH --mail-type=END,FAIL' + '\n') 
			f.write('#SBATCH --mail-user=' + self.email + '\n')
			f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
			f.write('#SBATCH -o ' + 'foldx_buildmodel' + '.%N.%j.log' + '\n')
			f.write('#SBATCH -e ' + 'foldx_buildmodel' + '.%N.%j.log' + '\n')
			for foldx_buildmodel_mutation in self.additional_foldx_buildmodel_mutations:
				dirname = osp.join(self.foldx_buildmodel_additional_dir, foldx_buildmodel_mutation)
				f.write('./foldx_20231231 -f ' + osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg') + '\n')

	def get_additional_foldx_buildmodel_mutations(self):
		self.split_foldx_mutations()
		self.create_foldx_buildmodel_pdb_files()
		self.create_foldx_buildmodel_config_files()
		self.create_foldx_buildmodel_slurm_job('1')

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	parser.add_argument('-a', '--account') # def-*
	parser.add_argument('-f', '--foldx_executable_name') # foldx_20231231
	args = parser.parse_args()

	print('Getting additional FoldX BuildModel mutations...')
	c = AdditionalMutations(args.script_dir, args.scratch_dir, args.account, args.foldx_executable_name)
	c.get_additional_foldx_buildmodel_mutations()



if __name__=='__main__':
	main()