'''
Reorganizes and renames columns in *_mutation_edgotypes_quasi_null_wildtype.tsv
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import get_missense_info, write_missense_info_to_file, check_create_dir


class FinalEdgotypeFile:
	def __init__(self, names, data_dir, final_data_dir):
		self.names = names
		self.data_dir = data_dir
		self.final_data_dir = final_data_dir
		check_create_dir(final_data_dir)

	def rename_organize_columns_all(self):
		for name in self.names:
			print('Creating final mutation location/edgotype file for:', name)
			header_dict, missense_info = get_missense_info(osp.join(self.data_dir, name + '_mutation_edgotypes_quasi_null_wildtype.tsv'))
			final_missense_info_file = osp.join(self.final_data_dir, name + '_location_edgotype.tsv')
			self.rename_organize_columns(name, missense_info, header_dict, final_missense_info_file)
			# return

	def rename_organize_columns(self, name, missense_info, header_dict, final_missense_info_file):
		mutation_location_nums = [0, 0, 0] # num IR, num exposed-NIR, num buried-NIR
		mutation_edgotype_nums = [0, 0, 0] # num edgetic, num quasi-wildtype, num quasi-null
		final_missense_info = [] # num edgetic, num quasi-wildtype, num quasi-null
		# new header
		unchanged_header_columns = missense_info[0][:header_dict['on_interfacial_res']]
		new_header = unchanged_header_columns + ['on_interfacial_res', 'modeller_foldx_mutations', 'foldx_pssm_binding_ddg', 'PPI_perturbation',  'foldx_buildmodel_mutations',
										'foldx_buildmodel_folding_ddg', 'relative_solvent_accessibility', 'mutation_location', 'mutation_edgotype']
		# print(new_header)
		final_missense_info.append(new_header)
		# get mutation info
		for line in missense_info[1:]: 
			final_columns = []
			final_columns.extend(line[:header_dict['ir_mutation']])
			# only ir_mutation --> mutation_location & edgotype --> mutation_edgotype are moved to the end
			# and quasi_null_wildtype removed
			final_columns.extend(line[header_dict['modeller_foldx_mutations']:header_dict['edgotype']]) # 'modeller_foldx_mutations', 'foldx_pssm_binding_ddg', 'perturbations_physics'
			final_columns.extend(line[header_dict['foldx_buildmodel_mutations']:header_dict['quasi_null_wildtype']]) # 'foldx_buildmodel_mutations', 'foldx_buildmodel_folding_ddg', 'relative_solvent_accessibility'
			# get mutation_location column
			mutation_location = ''
			ir_mutation = line[header_dict['ir_mutation']]
			if ir_mutation == '1':
				mutation_location = 'IR'
				mutation_location_nums[0] += 1
			else: # exposed-NIR or buried-NIR, look at RSA
				# print(ir_mutation)
				rsa = float(line[header_dict['relative_solvent_accessibility']])
				if rsa > 0.25: # exposed mutation
					mutation_location = 'exposed-NIR'
					mutation_location_nums[1] += 1
				else: # buried mutation
					mutation_location = 'buried-NIR'
					mutation_location_nums[2] += 1
			final_columns.append(mutation_location)
			# get mutation_edgotype column
			mutation_edgotype = ''
			edgotype = line[header_dict['edgotype']]
			if edgotype == 'edgetic':
				mutation_edgotype = 'edgetic'
				mutation_edgotype_nums[0] += 1
			else: # quasi-wildtype or quasi-null
				quasi_type = line[header_dict['quasi_null_wildtype']]
				if quasi_type != 'quasi-null' and quasi_type != 'quasi-wildtype': # sanity check
					print('Something went wrong!', quasi_type)
				else:
					mutation_edgotype = quasi_type
				if quasi_type == 'quasi-wildtype':
					mutation_edgotype_nums[1] += 1
				else:
					mutation_edgotype_nums[2] += 1
			final_columns.append(mutation_edgotype)
			# add to final_missense_info
			if len(final_columns) != len(new_header):
				print('Oh no!! Incorrect number of columns!')
			else:
				final_missense_info.append(final_columns)
		# write to file
		write_missense_info_to_file(final_missense_info, final_missense_info_file)
		# check nums
		# print(mutation_location_nums)
		# print(mutation_edgotype_nums)


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed', 'edgotypes')
	final_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'edgotypes_final')


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
			else: # data == 'dbsnp'
				names.append('_'.join([interactome, data, 'missense', 'mutations']))


	# for name in names:
	# 	print(name)

	f = FinalEdgotypeFile(names, data_dir, final_data_dir)
	f.rename_organize_columns_all()

if __name__=='__main__':
	main()