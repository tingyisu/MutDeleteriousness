'''
Process .pdb files outputted by MODELLER for input into FoldX PSSM
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''


import os
import os.path as osp
import argparse
import pickle
from Bio.PDB import PDBIO, PDBParser
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump, get_missense_info, write_missense_info_to_file

# IDEA:
# combine the two .pdb files from MODELLER to create a 'complex' for input into FoldX
# this 'complex' only contains 2 chains
# protein with mutation will always be chain A, while interactor protein will be chain B
# then creates jobs files to run FoldX PSSM (to find IR mutation induced change in binding free energy)

# NOTE: no email will be sent when each SLURM job is completed as there are too many jobs and the jobs will be completed in the allocated time


class CreatePDB:
	def __init__(self, script_dir, scratch_dir, account, foldx_executable_name):
		self.account = account
		self.foldx_executable_name = foldx_executable_name
		self.data_dir = osp.join(script_dir, 'files')
		self.results_dir = osp.join(script_dir, 'results')
		self.foldx_pssm_dir = osp.join(scratch_dir, 'foldx_pssm_all')
		self.foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		check_create_dir(self.foldx_pssm_dir)
		check_create_dir(self.foldx_buildmodel_dir)
		self.uniprot_alignment_dict = pickle_load(osp.join(self.data_dir, 'uniprot_alignment_dict.pickle')) # key = (uniprot, pdb_chain), value = dict of modeller residue numbering (starting from one), and original uniprot residue numbering key-value pairs
		self.all_foldx_pssm_mutations = pickle_load(osp.join(self.data_dir,'all_foldx_pssm_mutations.pickle'))
		self.all_foldx_buildmodel_mutations = pickle_load(osp.join(self.data_dir,'all_foldx_buildmodel_mutations.pickle'))
		self.num_mutations_per_job = 2000
		self.foldx_pssm_num_splits = int(len(self.all_foldx_pssm_mutations)/2000) + 1
		self.foldx_buildmodel_num_splits = int(len(self.all_foldx_buildmodel_mutations)/2000) + 1

	def create_parts(self, num_splits, dirname, mutations, name):
		total_num = 0
		total = []
		for i in range(num_splits):
			check_create_dir(osp.join(dirname, 'split_' + str(i)))
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
		if not check_exist(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_0.pickle')):
			self.create_parts(self.foldx_pssm_num_splits, self.foldx_pssm_dir, self.all_foldx_pssm_mutations, 'foldx_pssm_mutations')
		if not check_exist(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_0.pickle')):
			self.create_parts(self.foldx_buildmodel_num_splits, self.foldx_buildmodel_dir, self.all_foldx_buildmodel_mutations, 'foldx_buildmodel_mutations')
				
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

	def create_foldx_pssm_pdb_files(self):
		parser = PDBParser(QUIET=True)
		io = PDBIO()
		for i in range(self.foldx_pssm_num_splits):
			print('Creating PDB files for FoldX PSSM mutations part:', i)
			foldx_pssm_mutations_for_part = pickle_load(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_' + str(i) + '.pickle'))
			num_mutations = 1
			for foldx_pssm_mutation in foldx_pssm_mutations_for_part:
				model1, model2, mutation = foldx_pssm_mutation.split('-')
				model1_pdb_file = osp.join(self.results_dir, model1 + '.B99990001.pdb')
				model2_pdb_file = osp.join(self.results_dir, model2 + '.B99990001.pdb')
				dirname = osp.join(self.foldx_pssm_dir, 'split_' + str(i), foldx_pssm_mutation)
				# create new directory for each foldx_pssm_mutation
				check_create_dir(dirname)
				combined_pdb_file = osp.join(dirname, foldx_pssm_mutation + '.pdb')
				
				# get model1 and model2 structures
				structure1 = parser.get_structure(model1, model1_pdb_file)
				structure2 = parser.get_structure(model2, model2_pdb_file)

				# rename residue numbers in both model1 & model2 to uniprot numbering
				# get residue list for model1
				residue1_id_list = []
				for residue in structure1[0]['A']:
					residue1_id_list.append(residue.id)

				# get residue id dict for model1
				# need to do this in order to not get value error with overlapping residue numbers (https://github.com/biopython/biopython/issues/1551)
				residue1_id_dict = self.get_residue_id_dict(model1, residue1_id_list)
				for residue_id in residue1_id_list[::-1]: # change residue number in reverse to avoid overlapping residue numbers
					structure1[0]['A'][residue_id].id = residue1_id_dict[residue_id]

				# get residue list for model2 while renaming chain 'A' to 'B'
				residue2_id_list = []
				structure2[0]['A'].id = 'B' # rename chain
				for residue in structure2[0]['B']:
					residue2_id_list.append(residue.id)

				# get residue id dict for model2
				residue2_id_dict = self.get_residue_id_dict(model2, residue2_id_list)
				for residue_id in residue2_id_list[::-1]: # change residue number in reverse to avoid overlapping residue numbers
					structure2[0]['B'][residue_id].id = residue2_id_dict[residue_id]

				structure2[0]['B'].detach_parent()
				structure1[0].add(structure2[0]['B']) # combine chain B in model2 to model1

				io.set_structure(structure1)
				io.save(combined_pdb_file)
					
				print('Done:', str(num_mutations) + '/' + str(len(foldx_pssm_mutations_for_part)))
				num_mutations += 1

	def create_foldx_pssm_config_files(self):
		print('Creating FoldX PSSM configuration files...')
		for i in range(self.foldx_pssm_num_splits):
			print('Creating FoldX PSSM configuration files for part:', i)
			foldx_pssm_mutations_for_part = pickle_load(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_' + str(i) + '.pickle'))
			for foldx_pssm_mutation in foldx_pssm_mutations_for_part:
				dirname = osp.join(self.foldx_pssm_dir, 'split_' + str(i), foldx_pssm_mutation)
				_, _, foldx_mutation = foldx_pssm_mutation.split('-')
				mut_res = foldx_mutation[-1]
				position = foldx_mutation[:-1] + 'a'
				# write Pssm config file
				with open(osp.join(dirname, foldx_pssm_mutation + '_pssm_config.cfg'), 'w') as f:
					f.write('command=Pssm' + '\n')
					f.write('analyseComplexChains=A,B' + '\n')
					f.write('aminoacids=' + mut_res + '\n')
					f.write('positions=' + position + '\n')
					f.write('pdb-dir=' + dirname + '\n')
					f.write('output-dir=' + dirname + '\n')
					f.write('pdb=' + foldx_pssm_mutation + '.pdb' + '\n')
					# the rest of the parameters are all default

	def create_foldx_pssm_slurm_job(self, num_days):
		for i in range(self.foldx_pssm_num_splits):
			print('Creating FoldX PSSM SLURM job for part:', i)
			with open(osp.join(self.foldx_pssm_dir, 'split_' + str(i), 'foldx_pssm_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -N 1' + '\n')
				f.write('#SBATCH -n 5' + '\n')
				f.write('#SBATCH --mem 20G' + '\n')
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
				f.write('#SBATCH -o ' + 'foldx_pssm' + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e ' + 'foldx_pssm' + '.%N.%j.log' + '\n')
				foldx_pssm_mutations_for_part = pickle_load(osp.join(self.foldx_pssm_dir, 'foldx_pssm_mutations_' + str(i) + '.pickle'))
				for foldx_pssm_mutation in foldx_pssm_mutations_for_part:
					dirname = osp.join(self.foldx_pssm_dir, 'split_' + str(i), foldx_pssm_mutation)
					f.write('./' + self.foldx_executable_name + ' -f ' + osp.join(dirname, foldx_pssm_mutation + '_pssm_config.cfg') + '\n')
					# remove all intermediate/uneeded output files to conserve space (# of files) in scratch directory (maximum # of files = 1000k)
					f.write('rm ' + osp.join(dirname, '*.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'PSSM_Clash_'+ foldx_pssm_mutation + '.txt') + '\n')
					f.write('rm ' + osp.join(dirname, foldx_pssm_mutation + '_1.pdb') + '\n')
					f.write('rm ' + osp.join(dirname, 'WT_' + foldx_pssm_mutation + '_1.pdb') + '\n')
					f.write('rm ' + osp.join(dirname, 'individual_list_0_PSSM.txt') + '\n')

	def create_foldx_buildmodel_pdb_files(self):
		parser = PDBParser(QUIET=True)
		io = PDBIO()
		for i in range(self.foldx_buildmodel_num_splits):
			print('Creating PDB files for FoldX BuildModel mutations part:', i)
			foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_' + str(i) + '.pickle'))
			num_mutations = 1
			for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
				model, mutation = foldx_buildmodel_mutation.split('-')
				model_pdb_file = osp.join(self.results_dir, model + '.B99990001.pdb')
				dirname = osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), foldx_buildmodel_mutation)
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
		for i in range(self.foldx_buildmodel_num_splits):
			print('Creating FoldX BuildModel configuration files for part:', i)
			foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_' + str(i) + '.pickle'))
			for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
				dirname = osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), foldx_buildmodel_mutation)
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
		for i in range(self.foldx_buildmodel_num_splits):
			print('Creating FoldX BuildModel SLURM job for part:', i)
			with open(osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), 'foldx_buildmodel_' + str(i) + '.slurm'), 'w') as f:
				f.write('#!/bin/bash' + '\n')
				f.write('#SBATCH -N 1' + '\n')
				f.write('#SBATCH -n 5' + '\n')
				f.write('#SBATCH --mem 20G' + '\n')
				f.write('#SBATCH --account=' + self.account + '\n')
				f.write('#SBATCH -t ' + str(num_days) + '-00:00' + '\n')
				f.write('#SBATCH -o ' + 'foldx_buildmodel' + '.%N.%j.log' + '\n')
				f.write('#SBATCH -e ' + 'foldx_buildmodel' + '.%N.%j.log' + '\n')
				foldx_buildmodel_mutations_for_part = pickle_load(osp.join(self.foldx_buildmodel_dir, 'foldx_buildmodel_mutations_' + str(i) + '.pickle'))
				for foldx_buildmodel_mutation in foldx_buildmodel_mutations_for_part:
					dirname = osp.join(self.foldx_buildmodel_dir, 'split_' + str(i), foldx_buildmodel_mutation)
					f.write('./' + self.foldx_executable_name + ' -f ' + osp.join(dirname, foldx_buildmodel_mutation + '_buildmodel_config.cfg') + '\n')
					# remove all intermediate/uneeded output files to conserve space (# of files) in scratch directory (maximum # of files = 1000k)
					f.write('rm ' + osp.join(dirname, 'Average_' + foldx_buildmodel_mutation + '.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'PdbList_' + foldx_buildmodel_mutation + '.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, foldx_buildmodel_mutation + '_1.pdb') + '\n')
					f.write('rm ' + osp.join(dirname, 'Raw_' + foldx_buildmodel_mutation + '.fxout') + '\n')
					f.write('rm ' + osp.join(dirname, 'WT_' + foldx_buildmodel_mutation + '_1.pdb') + '\n')


	# FoldX PSSM
	def prepare_modeller_foldx_pssm_pdb_files(self):
		self.create_foldx_pssm_pdb_files()
		self.create_foldx_pssm_config_files()
		self.create_foldx_pssm_slurm_job(1)
	
	# FoldX BuildModel
	def prepare_modeller_foldx_buildmodel_pdb_files(self):
		self.create_foldx_buildmodel_pdb_files()
		self.create_foldx_buildmodel_config_files()
		self.create_foldx_buildmodel_slurm_job(1)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	parser.add_argument('-a', '--account') # def-* (compute canada project allocation/group name)
	parser.add_argument('-f', '--foldx_executable_name') # foldx_20231231
	args = parser.parse_args()

	print('Preparing modeller PDB files for Calculating Binding DDG (FoldX PSSM) and Folding DDG (FoldX BuildModel)')
	c = CreatePDB(args.script_dir, args.scratch_dir, args.account, args.foldx_executable_name)
	c.split_foldx_mutations()
	c.prepare_modeller_foldx_pssm_pdb_files()
	c.prepare_modeller_foldx_buildmodel_pdb_files()


if __name__=='__main__':
	main()