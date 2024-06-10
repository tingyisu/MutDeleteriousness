'''
Selects the PDB structural template for each protein-protein interaction
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import pickle
import time
import requests
from simple_tools import pickle_dump, pickle_load, check_exist

'''
* thresholds bitscore & resolution 
* Resolution for PDB structure must be >0 and <= 3.5 AND bitscore >= n (choose n, here set to 50)
* Once a PDB structural template meets the above, then select the one with the largest avg(bitscore)
* If more than one PDB structural temphate has the max average bitscore, then select the one with the longest combined interfacial residues
* the mmCIF files of the PDB structural templates are used along with their corresponding Uniprot protein sequences in MODELLER to build a 3D structural model from the Uniprot sequence alone
* here we don't need to worry about the residue discrepancies between Uniprot and PDB because MODELLER uses the Uniprot sequence in the resulting 3D model
'''

class SelectPDBTemplate:
	def __init__(self, name, data_orig_dir, data_processed_dir, resolution_threshold, bitscore_threshold):
		self.name = name
		self.resolution_threshold = resolution_threshold
		self.bitscore_threshold = bitscore_threshold
		self.resolutions_file = osp.join(data_orig_dir, 'resolu.idx')
		self.resolutions_dict_file = osp.join(data_processed_dir, 'interactome', 'pdb_resolutions_dict.pickle')
		# files to be read
		if name == 'intact':
			self.residues_file = osp.join(data_processed_dir, 'interactome', name + '_final_interfacial_residues_physical.tsv')
			self.blast_file = osp.join(data_processed_dir, 'interactome', name + '_blast_best_alignments_physical.tsv')
		else:
			self.residues_file = osp.join(data_processed_dir, 'interactome', name + '_final_interfacial_residues.tsv')
			self.blast_file = osp.join(data_processed_dir, 'interactome', name + '_blast_best_alignments.tsv')
		# files to be written to
		self.selected_interfacial_residues_file = osp.join(data_processed_dir, 'interactome', name + '_interfacial_residues_for_selected_template.tsv') # 
		self.selected_interfacial_residues_pickle_file = osp.join(data_processed_dir, 'interactome', name + '_interfacial_residues_for_selected_template.pickle')
		self.selected_pdb_structural_template_pickle_file = osp.join(data_processed_dir, 'interactome', name + '_selected_pdb_structural_template.pickle')
	
		self.individual_interfacial_residues = {} # dict of dict of interfacial residues for each PDB structural template (these interfacial residues are in Uniprot numbering)
		'''
		{(Protein1, Protein2):
			{
				'4awl_A_C': [residue1_list, residue2_list] 
			}

		}
		'''
		self.pdb_structures = [] # list of all pdb structures, used to find resolutions
		self.pdb_resolutions = {} # key = pdb structure, value = highest resolution (if have more than one resolution)
		self.blast_bitscore = {} # dict of dict of bitscore for each protein
		'''
		{protein:
			{
				'4awl_A': average bitscore
			}
		}
		'''
		self.bitscore_dict = {} # dict of dict of bitscores for protein pairs
		'''
		{(Protein1, Protein2):
			{
				'4awl_A_C': average bitscore
			}
		}
		'''
		self.interfacial_residues = {} # key = protein pair, value = list of interfacial residues (for selected pdb structural template, in Uniprot numbering)
		'''
		{(Protein1, Protein2): 
			[[123, 126, 159, 178], [55, 56, 59, 60]]
			...
		'''
		self.selected_pdb_templates = {} # key = protein pair, value = pdb structural template
		self.uniprot_proteins_pickle_file = osp.join(data_processed_dir, 'mutations', name + '_uniprot_proteins_to_map_mutations_to.pickle') # list of all UniProt proteins (to be used to map mutations to UniProt proteins)

	def convert_res_str_to_int(self, res):
		res_str = res.split(',')
		res_int = [int(res) for res in res_str]
		return res_int

	def update_individual_interfacial_residues_all_interaction_pairs(self, pdb, chain1, chain2, p1, p2, res1, res2):
		# update self.individual_interfacial_residues
		if (p2, p1) not in self.individual_interfacial_residues: # don't add reverse interaction
			if (p1, p2) not in self.individual_interfacial_residues:
				self.individual_interfacial_residues[(p1, p2)] = {}
			self.individual_interfacial_residues[(p1, p2)][pdb + '_' + chain1 + '_' + chain2] = [self.convert_res_str_to_int(res1), self.convert_res_str_to_int(res2)]

	# updates self.individual_interfacial_residues and self.pdb_structures
	def get_residues(self):
		print('-----Getting interfacial residues for each PPI structural model-----')
		with open(self.residues_file, 'r') as f:
			next(f)
			if self.name == 'hiunion':
				for line in f:
					items = line.strip().split('\t')
					p1, p2, pdb, chain1, chain2 = items[2:7]
					res1, res2 = items[7:9]
					self.update_individual_interfacial_residues_all_interaction_pairs(pdb, chain1, chain2, p1, p2, res1, res2)

			else: # self.name == 'intact'
				for line in f:
					items = line.strip().split('\t')
					p1, p2, pdb, chain1, chain2 = items[4:9]
					res1, res2 = items[9:11]
					self.update_individual_interfacial_residues_all_interaction_pairs(pdb, chain1, chain2, p1, p2, res1, res2)


	# get resolutions that meet the threshold (>0.0 AND <= self.resolution_threshold)
	# if a PDB structure has multiple resolutions, then select the smallest one in angstroms (that's not -1.0)
	def get_resolutions_threshold(self):
		# if not check_exist(self.resolutions_dict_file):
		with open(self.resolutions_file, 'r') as f_resolutions:
			for line in f_resolutions:
				if '\t' in line:
					items = line.strip().split('\t')
					# print(items)
					# break
					if len(items) == 3:
						pdb, _, resolution = items
						pdb = pdb.lower() # get lowercase pdb name
						if pdb != '':
							resolution = float(resolution) # convert resolution to float
							if resolution != -1.0 and resolution > 0.0 and resolution <= self.resolution_threshold: 
								# resolution doesn't exist if = -1.0
								# add all resolutions to self.pdb_resolutions_list
								if pdb in self.pdb_resolutions:
									prev_resolution = self.pdb_resolutions[pdb]
									if resolution < prev_resolution: # update with the smallest one
										self.pdb_resolutions[pdb] = resolution
								else:
									self.pdb_resolutions[pdb] = resolution
						else:
							print('Something went wrong with the pdb structure and its resolution...', pdb, resolution)

			pickle_dump(self.pdb_resolutions, self.resolutions_dict_file)

	# updates self.bitscore_dict with bitscores that reach a certain threshold
	def get_bitscore_dict_threshold(self):
		# get individual bitscores
		with open(self.blast_file, 'r') as f_blast:
			for line in f_blast:
				items = line.strip().split('\t')
				p, pdb, bitscore = items[0], items[1], float(items[11]) # here pdb = the structural template
				pdb_structure = pdb.split('_')[0]
				if pdb_structure in self.pdb_resolutions: # only take alignments of PDB structures that meet the resolution threshold
					if bitscore >= self.bitscore_threshold:
						if p not in self.blast_bitscore:
							self.blast_bitscore[p] = {}
						self.blast_bitscore[p][pdb] = bitscore
		# get average bitscores for each PPI structural template in a pair of proteins
		for pair in self.individual_interfacial_residues:
			for pdb_structural_template in self.individual_interfacial_residues[pair]:
				p1, p2 = pair
				pdb, chain1, chain2 = pdb_structural_template.split('_')
				pdb1, pdb2 = pdb + '_' + chain1, pdb + '_' + chain2
				if p1 in self.blast_bitscore and p2 in self.blast_bitscore:
					if pdb1 in self.blast_bitscore[p1] and pdb2 in self.blast_bitscore[p2]:
						bitscore1, bitscore2 = self.blast_bitscore[p1][pdb1], self.blast_bitscore[p2][pdb2]
						if pair not in self.bitscore_dict:
							self.bitscore_dict[pair] = {}
						self.bitscore_dict[pair][pdb_structural_template] = (bitscore1 + bitscore2)/2

	# get the pdb structural template with the longest total number of interfacial residues (would ideally have better coverage of edgetic mutations)
	# if have more than one, picks the first pdb structural template encountered
	def find_pdb_template_with_longest_residues(self, pair, pdb_templates_max_bitscore):
		max_residues_length, max_pdb_template = 0, ''
		for pdb_template in pdb_templates_max_bitscore:
			res1, res2 = self.individual_interfacial_residues[pair][pdb_template]
			if len(res1) + len(res2) > max_residues_length:
				max_residues_length = len(res1) + len(res2)
				max_pdb_template = pdb_template
		return max_pdb_template


	# updates self.interfacial_residues with the residues list of the selected PDB structural template
	# writes the residues list of the selected PDB structural template to file
	def pick_pdb_structural_templates(self):
		num_multiple_max_bitscore = 0
		# for each protein pair, see if there are more than one structural templates with the largest average bitscore
		for pair in self.bitscore_dict:
			bitscores = list(self.bitscore_dict[pair].values())
			# if pair == ('P27797', 'P30101'):
			# 	print(bitscores)
			keys = list(self.bitscore_dict[pair].keys())
			max_bitscore = max(bitscores)
			# get pdb structural templates with max_bitscore
			pdb_templates_max_bitscore = []
			for i in range(len(bitscores)):
				if bitscores[i] == max_bitscore:
					pdb_templates_max_bitscore.append(keys[i])
			# pick template with largest avg bitscore
			# if have more than one pdb structural template with max_bitscore
			# then pick the one with the longest number of residues
			if len(pdb_templates_max_bitscore) > 1:
				max_pdb_template = self.find_pdb_template_with_longest_residues(pair, pdb_templates_max_bitscore)
				self.selected_pdb_templates[pair] = max_pdb_template
				self.interfacial_residues[pair] = self.individual_interfacial_residues[pair][max_pdb_template]
				num_multiple_max_bitscore += 1
			else:
				# print('length of pdb_templates_max_bitscore:', len(pdb_templates_max_bitscore))
				if len(pdb_templates_max_bitscore) != 1:
					print('ERROR! len(pdb_templates_max_bitscore) != 1:', pdb_templates_max_bitscore)
				best_pdb_template = pdb_templates_max_bitscore[0]
				self.selected_pdb_templates[pair] = best_pdb_template
				self.interfacial_residues[pair] = self.individual_interfacial_residues[pair][best_pdb_template]

		# pickle self.interfacial_residues self.selected_pdb_templates
		pickle_dump(self.interfacial_residues, self.selected_interfacial_residues_pickle_file)
		pickle_dump(self.selected_pdb_templates, self.selected_pdb_structural_template_pickle_file)

		print('Number of PPIs with multiple max-bitscore templates:', num_multiple_max_bitscore)


	# writes interfacial residues for selected pdb structural templates to file
	def write_best_pdb_structural_template(self):
		with open(self.residues_file, 'r') as f, open(self.selected_interfacial_residues_file, 'w') as f_write:
			f_write.write(f.readline()) # write header line
			unique_uniprot_pairs = [] # only keep unique uniprot pairs, don't want multiples of the same pairs but with diff gene/ensembl id
			if self.name == 'hiunion':		
				for line in f:
					items = line.strip().split('\t')
					p1, p2, pdb, chain1, chain2 = items[2:7]
					pdb_structural_template = '_'.join([pdb, chain1, chain2])
					if (p1, p2) in self.selected_pdb_templates and (p1, p2) not in unique_uniprot_pairs: # take only first Ensembl pair for each PPI
						if pdb_structural_template in self.selected_pdb_templates[(p1, p2)]:
							f_write.write(line)
							unique_uniprot_pairs.append((p1, p2))
			else: # self.name == 'intact'
				for line in f:
					items = line.strip().split('\t')
					p1, p2, pdb, chain1, chain2 = items[4:9]
					pdb_structural_template = '_'.join([pdb, chain1, chain2])
					if (p1, p2) in self.selected_pdb_templates and (p1, p2) not in unique_uniprot_pairs: # take only first gene pair for each PPI
						if pdb_structural_template in self.selected_pdb_templates[(p1, p2)]:
							f_write.write(line)
							unique_uniprot_pairs.append((p1, p2))
			print('Number of unique uniprot pairs:', len(unique_uniprot_pairs))
			# count number of unique proteins
			unique_proteins = []
			for (p1, p2) in unique_uniprot_pairs:
				unique_proteins.extend([p1, p2])
			print('Number of unique proteins:', len(set(unique_proteins)))
			pickle_dump(unique_proteins, self.uniprot_proteins_pickle_file)

	def get_pdb_structures(self): # compiles all the PDB structures from the selected PDB structural templates
		for pair in self.selected_pdb_templates:
			pdb_structural_template = self.selected_pdb_templates[pair]
			pdb = pdb_structural_template.split('_')[0]
			# append PDB structure
			self.pdb_structures.append(pdb)

		# get only unique pdb structures
		self.pdb_structures = list(set(self.pdb_structures))

		print('There are', len(self.pdb_structures), 'unique PDB structures in the selected templates')
		

	def select_best_pdb_structural_templates(self):
		self.get_residues()
		self.get_resolutions_threshold()
		self.get_bitscore_dict_threshold()
		self.pick_pdb_structural_templates()
		self.write_best_pdb_structural_template()
		self.get_pdb_structures()

def main():
	script_dir = osp.dirname(__file__)
	data_orig_dir = osp.join(script_dir, '..', 'data', 'original')
	data_processed_dir = osp.join(script_dir, '..', 'data', 'processed')
	resolution_threshold = 3.5
	bitscore_threshold = 50

	print('-----HI-Union-----')

	hiunion = SelectPDBTemplate('hiunion', data_orig_dir, data_processed_dir, resolution_threshold, bitscore_threshold)
	hiunion.select_best_pdb_structural_templates()

	print('-----IntAct-----')

	intact = SelectPDBTemplate('intact', data_orig_dir, data_processed_dir, resolution_threshold, bitscore_threshold)
	intact.select_best_pdb_structural_templates()

if __name__=='__main__':
	main()