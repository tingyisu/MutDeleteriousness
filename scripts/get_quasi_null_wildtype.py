'''
Determine quasi-null and quasi-wildtype mutations using FoldX BuildModel and RSA calculations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import pickle_load, check_create_dir, get_missense_info, write_missense_info_to_file, check_exist
import argparse
import shutil
import math

class QuasiNullWildtype:
	def __init__(self, names, script_dir, scratch_dir):
		self.names = names
		self.data_dir = osp.join(script_dir, 'files')
		self.all_foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		self.foldx_buildmodel_combined_dir = osp.join(script_dir, 'folding_ddg')
		self.rsa_dir = osp.join(script_dir, 'rsa')
		check_create_dir(self.rsa_dir)
		check_create_dir(self.foldx_buildmodel_combined_dir)

	# move all buildmodel files from split_* and additional_*directories under self.foldx_buildmodel_dir into another directory
	def combine_copy_foldx_buildmodel_rsa_results(self):
		print('Copying FoldX BuildModel and RSA results from scratch directory...')
		split_additional_folders = [f.path for f in os.scandir(self.all_foldx_buildmodel_dir) if f.is_dir()]
		print(split_additional_folders, len(split_additional_folders))
		for folder in split_additional_folders:
			print('Copying RSA files for:', folder)
			mutation_folders = [f.path for f in os.scandir(folder) if f.is_dir() and 'molecules' not in f.path]
			# iterate through each mutation fiolder
			for mutation_folder in mutation_folders:
				mutation_name = mutation_folder.split('/')[-1]
				buildmodel_file = osp.join(mutation_folder, 'Dif_' + mutation_name + '.fxout')
				rsa_file = osp.join(mutation_folder, mutation_name + '.rsa.csv')
				# copy FoldX BuildModel file
				if check_exist(buildmodel_file):
					shutil.copy(buildmodel_file, self.foldx_buildmodel_combined_dir)
				# copy RSA file
				if check_exist(rsa_file):
					shutil.copy(rsa_file, self.rsa_dir)


	def get_quasi_null_wildtype_all(self):
		for name in self.names:
			print('Finding quasi null/wildtype mutations for:', name)
			# read in edgotypes file
			header_dict, missense_info = get_missense_info(osp.join(self.data_dir, name + '_mutation_edgotypes.tsv'))

			# get quasi-null/wildtype mutations and write to file
			missense_info_edgotype_quasi_null_wildtype = osp.join(self.data_dir, name + '_mutation_edgotypes_quasi_null_wildtype.tsv')
			self.get_quasi_null_wildtype(name, missense_info, header_dict, missense_info_edgotype_quasi_null_wildtype)

	def get_quasi_null_wildtype(self, name, missense_info, header_dict, missense_info_edgotype_quasi_null_wildtype):
		missense_info[0].append('foldx_buildmodel_folding_ddg')
		header_dict['foldx_buildmodel_folding_ddg'] = len(header_dict)
		missense_info[0].append('relative_solvent_accessibility')
		header_dict['relative_solvent_accessibility'] = len(header_dict)
		missense_info[0].append('quasi_null_wildtype')
		header_dict['quasi_null_wildtype'] = len(header_dict)

		for i in range(1, len(missense_info)):
			# print(header_dict)
			edgotype = missense_info[i][header_dict['edgotype']] # switch to edgotype?
			if edgotype == 'non-edgetic':
				foldx_buildmodel_mutation = missense_info[i][header_dict['foldx_buildmodel_mutations']]
				# print('Looking at:', foldx_buildmodel_mutation)
				mutation = foldx_buildmodel_mutation.split('-')[-1]
				mut_pos = ''.join(c for c in mutation if c.isdigit())
				folding_ddg = ''
				rsa = ''
				# get folding_ddg
				with open(osp.join(self.foldx_buildmodel_combined_dir, 'Dif_' + foldx_buildmodel_mutation + '.fxout'), 'r') as f_folding_ddg:
					folding_ddg = f_folding_ddg.readlines()[-1].split('\t')[1]
					missense_info[i].append(folding_ddg) # append folding_ddg
					folding_ddg = float(folding_ddg)
				# get rsa
				with open(osp.join(self.rsa_dir, foldx_buildmodel_mutation + '.rsa.csv'), 'r') as f_rsa:
					# actual_wt_res = foldx_buildmodel_mutation[-1]
					actual_wt_res = foldx_buildmodel_mutation.split('-')[1][0] # first character
					rsa_header = f_rsa.readline().strip().split(',')
					rsa_header_dict = {}
					j = 0
					for item in rsa_header:
						rsa_header_dict[item] = j
						j += 1
					for line in f_rsa:
						items = line.strip().split(',')
						position, wt_res, curr_rsa = items[rsa_header_dict['pdb_position']], items[rsa_header_dict['pdb_aa']], items[rsa_header_dict['rsa']]
						if position == mut_pos:
							if wt_res != actual_wt_res:
								print('incorrect mutation residue in RSA file:', foldx_buildmodel_mutation, '; have:', wt_res, ', should be:', actual_wt_res)
							else:
								missense_info[i].append(curr_rsa) # append RSA value
								rsa = float(curr_rsa)
								break
				# quasi-wildtype (exposed residue, mutation doesn't make protein fall apart)
				# quasi-null (buried residue, mutation makes protein fall apart)
				if rsa > 0.25: # exposed residue
					missense_info[i].append('quasi-wildtype')
				else: # buried residue (rsa <= 0.25)
					if folding_ddg >= 2.0: # mutation de-stabilizes protein
						missense_info[i].append('quasi-null')
					else: # folding_ddg < 2.0 (mutation isn't enough to de-stabilize protein)
						missense_info[i].append('quasi-wildtype')
			else: # edgetic mutation
				missense_info[i].append('NA')
				missense_info[i].append('NA')
				missense_info[i].append('NA')
		write_missense_info_to_file(missense_info, missense_info_edgotype_quasi_null_wildtype)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	args = parser.parse_args()

	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	dbsnp_pathogenicities = ['', 'pathogenic']
	names = []

	for interactome in interactomes:
		for data in mutation_data:
			if data == 'dbsnp':
				for dbsnp_pathogenicity in dbsnp_pathogenicities:
					if dbsnp_pathogenicity == '':
						names.append('_'.join([interactome, data, 'missense', 'mutations']))
					else:
						names.append('_'.join([interactome, data, 'missense', 'mutations', 'pathogenic']))
			else: # data == 'clinvar'
				names.append('_'.join([interactome, data, 'missense', 'mutations']))


	for name in names:
		print(name)


	print('-----Finding quasi null, quasi wildtype based on RSA and folding DDG-----')
	q = QuasiNullWildtype(names, args.script_dir, args.scratch_dir)
	q.combine_copy_foldx_buildmodel_rsa_results() # copying files over will take a long time, so need a .slurm file
	q.get_quasi_null_wildtype_all()

	

if __name__=='__main__':
	main()


