'''
Use FoldX PSSM (mutation induced change in binding free energies to find edgetic mutations)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import pickle_load, check_create_dir, get_missense_info, write_missense_info_to_file, pickle_dump, check_exist
import argparse
import shutil

class BindingDDG:
	def __init__(self, names, script_dir, scratch_dir):
		self.names = names
		self.data_dir = osp.join(script_dir, 'files')
		self.blast_best_alignments =  pickle_load(osp.join(self.data_dir, 'all_blast_best_alignments.pickle'))
		self.all_modeller_templates = pickle_load(osp.join(self.data_dir, 'all_modeller_templates.pickle'))
		self.all_foldx_buildmodel_mutations = pickle_load(osp.join(self.data_dir, 'all_foldx_buildmodel_mutations.pickle'))
		self.foldx_pssm_dir = osp.join(scratch_dir, 'foldx_pssm_all')
		self.foldx_pssm_combined_dir = osp.join(script_dir, 'binding_ddg')
		check_create_dir(self.foldx_pssm_combined_dir)
		self.foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		self.additional_modeller_templates = []
		self.additional_modeller_templates_pickle_file = osp.join(self.foldx_buildmodel_dir, 'additional_modeller_templates.pickle')
		self.additional_foldx_buildmodel_mutations = []
		self.additional_foldx_buildmodel_mutations_pickle_file = osp.join(self.foldx_buildmodel_dir, 'additional_foldx_buildmodel_mutations.pickle')

	# move all PSSM files from split_* directories under self.foldx_pssm_dir into another directory
	def combine_copy_foldx_pssm(self):
		split_folders = [f.path for f in os.scandir(self.foldx_pssm_dir) if f.is_dir()]
		# print(split_folders)
		for folder in split_folders:
			print('Copying PSSM files for:', folder)
			mutation_folders = [f.path for f in os.scandir(folder) if f.is_dir() and 'molecules' not in f.path]
			# print(mutation_folders)
			for mutation_folder in mutation_folders:
				mutation_name = mutation_folder.split('/')[-1]
				pssm_file = osp.join(mutation_folder, 'PSSM_' + mutation_name + '.txt')
				# print(mutation_name)
				if check_exist(pssm_file): # so that don't get an error if FoldX PSSM didn't work...
					shutil.copy(pssm_file, self.foldx_pssm_combined_dir)
			# break
		# shutil.copy(src, dst)

	# takes list of possible templates for uniprot protein with non-edgetic mutation
	# returns best template with the lowest evalue (if more than one, picks the first one)
	def select_non_edgetic_template(self, non_edgetic_templates):
		evalue_list = []
		for (uniprot_protein, pdb_chain) in non_edgetic_templates:
			alignment = self.blast_best_alignments[(uniprot_protein, pdb_chain)]
			evalue = float(alignment[8])
			evalue_list.append(evalue)
		min_evalue = min(evalue_list)
		min_index = evalue_list.index(min_evalue)
		# print(evalue_list, min_evalue, min_index, non_edgetic_templates)
		return non_edgetic_templates[min_index][1] # return pdb_chain

	def get_edgotype_all(self):
		for name in self.names:
			header_dict, missense_info = get_missense_info(osp.join(self.data_dir, name + '_ir_with_modeller_foldx_mutations.tsv'))
			missense_info_edgotype = osp.join(self.data_dir, name + '_mutation_edgotypes.tsv')
			interactome_name = name.split('_')[0]
			selected_pdb_structural_template = pickle_load(osp.join(self.data_dir, interactome_name + '_selected_pdb_structural_template.pickle'))
			self.get_binding_ddg_edgotype(name, missense_info, header_dict, missense_info_edgotype, selected_pdb_structural_template)
		print('-----Additional items to run-----')
		self.additional_modeller_templates = list(set(self.additional_modeller_templates))
		self.additional_foldx_buildmodel_mutations = list(set(self.additional_foldx_buildmodel_mutations))
		print('Number of additional MODELLER templates to run:', len(self.additional_modeller_templates))
		print('Number of additional FoldX BuildModel mutations:', len(self.additional_foldx_buildmodel_mutations))
		pickle_dump(self.additional_modeller_templates, self.additional_modeller_templates_pickle_file)
		pickle_dump(self.additional_foldx_buildmodel_mutations, self.additional_foldx_buildmodel_mutations_pickle_file)

	def get_binding_ddg_edgotype(self, name, missense_info, header_dict, missense_info_edgotype, selected_pdb_structural_template):
		edgetic, non_edgetic = 0, 0
		edgetic_turned_non_edgetic = 0
		missense_info[0].append('foldx_pssm_binding_ddg')
		header_dict['foldx_pssm_binding_ddg'] = len(header_dict)
		missense_info[0].append('perturbations')
		header_dict['perturbations'] = len(header_dict)
		missense_info[0].append('edgotype')
		header_dict['edgotype'] = len(header_dict)
		missense_info[0].append('foldx_buildmodel_mutations')
		header_dict['foldx_buildmodel_mutations'] = len(header_dict)
		for i in range(1, len(missense_info)):
			ir_mutation = missense_info[i][header_dict['ir_mutation']]
			on_interfacial_res = missense_info[i][header_dict['on_interfacial_res']]
			if ir_mutation == '0': # not on interfacial residue
				modeller_foldx_mutations = missense_info[i][header_dict['modeller_foldx_mutations']]
				missense_info[i].append('NA')
				missense_info[i].append(on_interfacial_res)
				missense_info[i].append('non-edgetic')
				missense_info[i].append(modeller_foldx_mutations)
				non_edgetic += 1
			else: # on interfacial residue
				modeller_foldx_mutations = missense_info[i][header_dict['modeller_foldx_mutations']].split(',')
				binding_ddgs = []
				perturbations = []
				foldx_mutation = ''
				for modeller_foldx_mutation in modeller_foldx_mutations:
					if modeller_foldx_mutation != '-1':
						pssm_file = osp.join(self.foldx_pssm_combined_dir, 'PSSM_' + modeller_foldx_mutation + '.txt')
						# get foldx mutation
						if foldx_mutation == '':
							foldx_mutation = modeller_foldx_mutation.split('-')[-1]
						with open(pssm_file, 'r') as f:
							lines = f.readlines()
							# mut_res = lines[0].strip()
							_, ddg = lines[1].strip().split('\t')
							binding_ddgs.append(ddg)
							ddg = float(ddg)
							# print(mut_res, position, ddg)
							if ddg > 0.5: # if mutation is at least mildly destabilizing
								perturbations.append('1')
							else: # ddg <= 0.5
								perturbations.append('0')
					else:
						binding_ddgs.append('NA')
						perturbations.append('0')
				missense_info[i].append(','.join(binding_ddgs))
				missense_info[i].append(','.join(perturbations))
				if '1' in perturbations:
					missense_info[i].append('edgetic')
					missense_info[i].append('NA')
					edgetic += 1
				else:
					missense_info[i].append('non-edgetic')
					non_edgetic += 1
					edgetic_turned_non_edgetic += 1
					non_edgetic_templates = []
					uniprot_protein = missense_info[i][header_dict['uniprot_protein']]
					protein_partners = missense_info[i][header_dict['protein_partners']].split(',')
	
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

						# add non_edgetic templates
						if (uniprot_protein, pdb_chain) not in non_edgetic_templates:
							non_edgetic_templates.append((uniprot_protein, pdb_chain))

					# get unique templates
					if len(non_edgetic_templates) == 1:
						best_template = non_edgetic_templates[0][1]
					else:
						# find the template with the smallest evalue alignment (if more than one, takes first occurrence)
						best_template = self.select_non_edgetic_template(non_edgetic_templates)

					model = '_'.join([uniprot_protein, best_template])
					if model not in self.all_modeller_templates:
						self.additional_modeller_templates.append(model)
					foldx_buildmodel_mutation = '-'.join([model, foldx_mutation])
					# add foldx_build_model_mutation to self.missense_info
					missense_info[i].append(foldx_buildmodel_mutation)
					if foldx_buildmodel_mutation not in self.all_foldx_buildmodel_mutations:
						self.additional_foldx_buildmodel_mutations.append(foldx_buildmodel_mutation) # add to list, to pickle later
		
		print('-----' + name + '-----')
		print('Number of edgetic mutations:', edgetic)
		print('Number of non-edgetic mutations:', non_edgetic)
		print('Number of edgetic turned non-edgetic mutations:', edgetic_turned_non_edgetic)
		write_missense_info_to_file(missense_info, missense_info_edgotype)


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	args = parser.parse_args()


	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	names = []

	for interactome in interactomes:
		for data in mutation_data:
			names.append('_'.join([interactome, data, 'missense', 'mutations']))


	print('-----Finding edgotypes (edgetic/non-edgetic) using FoldX binding DDG-----')
	b = BindingDDG(names, args.script_dir, args.scratch_dir)
	b.combine_copy_foldx_pssm()
	b.get_edgotype_all()

	

if __name__=='__main__':
	main()


