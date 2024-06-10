'''
Combine BLASTP alignments from HI-Union and IntAct interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump, get_missense_info, write_missense_info_to_file

# get templates that need to be modelled with MODELLER
# these templates are either templates for disrupted interactions (in IR mutations) or from a specifically selected interaction 
# on a protein with a non-IR mutation

class ModellerTemplates:
	def __init__(self, data_dir, mutation_files_list):
		self.data_dir = data_dir
		self.mutation_files_list = mutation_files_list
		self.blast_best_alignments =  pickle_load(osp.join(self.data_dir, 'interactome', 'all_blast_best_alignments.pickle'))
		self.modeller_templates_pickle_file = osp.join(self.data_dir, 'edgotypes', 'all_modeller_templates.pickle')
		self.foldx_pssm_mutations_pickle_file = osp.join(self.data_dir, 'edgotypes', 'all_foldx_pssm_mutations.pickle')
		self.foldx_buildmodel_mutations_pickle_file = osp.join(self.data_dir, 'edgotypes', 'all_foldx_buildmodel_mutations.pickle')
		self.pdb_structures_pickle_file = osp.join(self.data_dir, 'edgotypes', 'pdb_structures_to_download.pickle')
		self.modeller_templates = [] # all templates that need to be modelled with MODELLER
		self.foldx_pssm_mutations = [] # all mutations on which to run FoldX PSSM
		self.foldx_buildmodel_mutations = [] # all mutations on which to run FoldX BuildModel
		self.pdb_structures = []
		# for mutation positions that are over 10000 (cannot be accomodated by .pdb files, so need to take %10000)
		self.pssm_over_9999 = [] # list of tuples (foldx pssm mutation, original mutation position prior to taking %)
		self.buildmodel_over_9999 = [] # list of tuples (foldx pssm mutation, original mutation position prior to taking %)

	def get_all_templates_mutations(self):
		for ir_mutations_file in self.mutation_files_list:
			print('Getting modeller templates and foldx mutations for:', ir_mutations_file)

			# get necessary info
			interactome = ir_mutations_file.split('_')[0]
			selected_pdb_structural_template = pickle_load(osp.join(self.data_dir, 'interactome', interactome + '_selected_pdb_structural_template.pickle'))
			header_dict, missense_info = get_missense_info(osp.join(self.data_dir, 'edgotypes', ir_mutations_file))
			if len(header_dict) != len(missense_info[1]): # sanity check
				print('Header & entries have different lengths!', name, len(header_dict), len(missense_info[1]))

			# file to write MODELLER & FoldX templates/mutations to 
			ir_modeller_foldx_mutations_file = osp.join(self.data_dir, 'edgotypes', ir_mutations_file.split('.')[0] + '_with_modeller_foldx_mutations.tsv')

			# updates self.modeller_templates, self.foldx_pssm_mutations, and self.foldx_buildmodel_mutations for each ir_mutations_file
			self.get_modeller_templates_foldx_mutations(header_dict, missense_info, selected_pdb_structural_template)
			# write missense info with included modeller foldx mutations to file
			write_missense_info_to_file(missense_info, ir_modeller_foldx_mutations_file)


		# take only unique entries in each
		self.modeller_templates = list(set(self.modeller_templates))
		self.foldx_pssm_mutations = list(set(self.foldx_pssm_mutations))
		self.foldx_buildmodel_mutations = list(set(self.foldx_buildmodel_mutations))
		self.pssm_over_9999 = list(set(self.pssm_over_9999))
		self.buildmodel_over_9999 = list(set(self.buildmodel_over_9999))
		print('Total number of modeller_templates:', len(self.modeller_templates))
		print('Total number of FoldX PSSM mutations:', len(self.foldx_pssm_mutations))
		print('Total number of foldx_buildmodel_mutations:', len(self.foldx_buildmodel_mutations))
		print('PSSM mutations with positions over 9999:', self.pssm_over_9999)
		print('BuildModel mutations with positions over 9999:', self.buildmodel_over_9999)
		# get list of pdb structures needed for MODELLER templates
		for template in self.modeller_templates:
			uniprot_protein, pdb, chain = template.split('_')
			self.pdb_structures.append(pdb)
		self.pdb_structures = list(set(self.pdb_structures))
	
		print('Total number of PDB structures to download:', len(self.pdb_structures))
		# pickle into corresponding files
		pickle_dump(self.modeller_templates, self.modeller_templates_pickle_file)
		pickle_dump(self.foldx_pssm_mutations, self.foldx_pssm_mutations_pickle_file)
		pickle_dump(self.foldx_buildmodel_mutations, self.foldx_buildmodel_mutations_pickle_file)
		pickle_dump(self.pdb_structures, self.pdb_structures_pickle_file)


	# takes list of possible templates for uniprot protein with non-ir mutation
	# returns best template with the lowest evalue (if more than one, picks the first one)
	def select_non_ir_template(self, non_ir_templates):
		evalue_list = []
		for (uniprot_protein, pdb_chain) in non_ir_templates:
			alignment = self.blast_best_alignments[(uniprot_protein, pdb_chain)]
			evalue = float(alignment[8])
			evalue_list.append(evalue)
		min_evalue = min(evalue_list)
		min_index = evalue_list.index(min_evalue)
		# print(evalue_list, min_evalue, min_index, non_ir_templates)
		return non_ir_templates[min_index][1] # return pdb_chain

	# mutation is ir if it maps onto an interfacial residue of at least one interaction
	def get_modeller_templates_foldx_mutations(self, header_dict, missense_info, selected_pdb_structural_template):
		# print(header_dict)
		# print('-----Compiling templates to be modelled using MODELLERs-----')
		missense_info[0].append('modeller_foldx_mutations')
		header_dict['modeller_foldx_mutations'] = len(header_dict)
		for i in range(1, len(missense_info)):
			# print(len(missense_info[i]), len(header_dict))
			over_9999 = False
			uniprot_protein = missense_info[i][header_dict['uniprot_protein']]
			ir_mutation = missense_info[i][header_dict['ir_mutation']]
			protein_partners = missense_info[i][header_dict['protein_partners']].split(',')
			on_ir = missense_info[i][header_dict['on_interfacial_res']].split(',')
			ref_res = missense_info[i][header_dict['ref_res']]
			alt_res = missense_info[i][header_dict['alt_res']]
			prot_res_pos_in_uniprot_protein = missense_info[i][header_dict['prot_res_pos_in_uniprot_protein']]
			# get foldx_mutation, take % if mutation position > 9999
			new_prot_res_pos_in_uniprot_protein = -1
			if int(prot_res_pos_in_uniprot_protein) > 9999:
				new_prot_res_pos_in_uniprot_protein = str(int(prot_res_pos_in_uniprot_protein) % 10000)
				over_9999 = True
			else:
				new_prot_res_pos_in_uniprot_protein = prot_res_pos_in_uniprot_protein
			foldx_mutation = ref_res + 'A' + new_prot_res_pos_in_uniprot_protein + alt_res 

			if ir_mutation == '1': # is an IR mutation
				foldx_pssm_mutations_list = [] # for each mutation (line in missense file)
				for j in range(len(protein_partners)):
					
					if on_ir[j] == '1':
						model1, model2 = '', ''
						if (uniprot_protein, protein_partners[j]) in selected_pdb_structural_template:
							pdb, chain1, chain2 = selected_pdb_structural_template[(uniprot_protein, protein_partners[j])].split('_')
							model1, model2 = '_'.join([uniprot_protein, pdb, chain1]), '_'.join([protein_partners[j], pdb, chain2])
						elif (protein_partners[j], uniprot_protein) in selected_pdb_structural_template:
							pdb, chain2, chain1 = selected_pdb_structural_template[(protein_partners[j], uniprot_protein)].split('_')
							model1, model2 = '_'.join([uniprot_protein, pdb, chain1]), '_'.join([protein_partners[j], pdb, chain2])
						else:
							print('ERROR!!! interaction pair is not in self.selected_templates_dict!', uniprot_protein, protein_partners[j])
						# .pdb formats cannot handle cases where prot_res_pos_in_uniprot_protein > 9999, so do the following
						# add modeller templates
						self.modeller_templates.extend([model1, model2])
						foldx_pssm_mutation = '-'.join([model1, model2, foldx_mutation])
						if over_9999:
							self.pssm_over_9999.append((foldx_pssm_mutation, prot_res_pos_in_uniprot_protein))
						foldx_pssm_mutations_list.append(foldx_pssm_mutation)
						self.foldx_pssm_mutations.append(foldx_pssm_mutation) # add to list, to pickle later
					
					else:
						foldx_pssm_mutations_list.append('-1')
				
				# add foldx_pssm_mutations_list to self.missense_info
				missense_info[i].append(','.join(foldx_pssm_mutations_list))

			else: # non-ir
				non_ir_templates = []
				for j in range(len(protein_partners)):

					pdb_chain = ''
					if (uniprot_protein, protein_partners[j]) in selected_pdb_structural_template:
						pdb, chain1, _ = selected_pdb_structural_template[(uniprot_protein, protein_partners[j])].split('_')
						pdb_chain = '_'.join([pdb, chain1])
					elif (protein_partners[j], uniprot_protein) in selected_pdb_structural_template:
						pdb, _, chain1 = selected_pdb_structural_template[(protein_partners[j], uniprot_protein)].split('_')
						pdb_chain = '_'.join([pdb, chain1])
					else:
						print('ERROR!!! interaction pair is not in self.selected_templates_dict!', uniprot_protein, protein_partners[j])

					# add non_ir templates
					if (uniprot_protein, pdb_chain) not in non_ir_templates:
						non_ir_templates.append((uniprot_protein, pdb_chain))

				# get unique templates
				if len(non_ir_templates) == 1:
					best_template = non_ir_templates[0][1]
				else:
					# find the template with the smallest evalue alignment (if more than one, takes first occurrence)
					best_template = self.select_non_ir_template(non_ir_templates)
				model = '_'.join([uniprot_protein, best_template])
				self.modeller_templates.append(model)
				foldx_buildmodel_mutation = '-'.join([model, foldx_mutation])
				if over_9999:
					self.buildmodel_over_9999.append((foldx_buildmodel_mutation, prot_res_pos_in_uniprot_protein))
				# add foldx_build_model_mutation to self.missense_info
				missense_info[i].append(foldx_buildmodel_mutation)
				self.foldx_buildmodel_mutations.append(foldx_buildmodel_mutation) # add to list, to pickle later
		


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	dbsnp_pathogenicities = ['', 'pathogenic']
	mutation_files_list = []

	for interactome in interactomes:
		for data in mutation_data:
			if data == 'dbsnp':
				for dbsnp_pathogenicity in dbsnp_pathogenicities:
					if dbsnp_pathogenicity == '':
						mutation_files_list.append('_'.join([interactome, data, 'missense', 'mutations', 'ir.tsv']))
					else:
						mutation_files_list.append('_'.join([interactome, data, 'missense', 'mutations', 'pathogenic', 'ir.tsv']))
			else: # data == 'clinvar'
				mutation_files_list.append('_'.join([interactome, data, 'missense', 'mutations', 'ir.tsv']))

	templates = ModellerTemplates(data_dir, mutation_files_list)
	templates.get_all_templates_mutations()

if __name__=='__main__':
	main()