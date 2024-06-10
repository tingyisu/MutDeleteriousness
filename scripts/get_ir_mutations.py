'''
Finds dbSNP and ClinVar missense mutations that lie on interfacial residues (IR mutations)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import os.path as osp
import os
from simple_tools import pickle_load, get_missense_info, write_missense_info_to_file

# removes mutations that do not have protein partners or that don't have any leftoever interactions in which the mutation position lies within the structural template

class IRMutation:
	def __init__(self, data_dir, final_data_dir, mutation_files_list):
		self.data_dir = data_dir
		self.final_data_dir = final_data_dir
		self.mutation_files_list = mutation_files_list
	

	# gets interactors and initializes missense_info_with_protein_partners (are the missense mutations on proteins that exist in the list of interactions)
	# only keeps mutations with at least one existing protein partner (interaction)
	def get_interactors(self, header_dict, missense_info, interfacial_residues):
		interactors = {}
		missense_info_with_protein_partners = []

		missense_info_with_protein_partners.append(missense_info[0])
		for i in range(1, len(missense_info)):
			uniprot_protein = missense_info[i][header_dict['uniprot_protein']]
			if uniprot_protein not in interactors:
				# get interacting partners
				partners = []
				for (p1, p2) in interfacial_residues:
					if uniprot_protein == p1:
						partners.append(p2)
					elif uniprot_protein == p2:
						partners.append(p1)
				interactors[uniprot_protein] = partners
				# add if has interacting partners
				if partners != []:
					missense_info_with_protein_partners.append(missense_info[i])
			else:
				# add if has interacting partners
				if interactors[uniprot_protein] != []:
					missense_info_with_protein_partners.append(missense_info[i])

		return interactors, missense_info_with_protein_partners

	# mutation is an IR mutation if it maps onto an interfacial residue of at least one interaction
	def find_ir_mutations(self, mutation_file, ir_mutations_file):
		header_dict, missense_info = get_missense_info(osp.join(self.final_data_dir, mutation_file))

		# get necessary info
		best_blast_alignment = {}
		interfacial_residues = {}
		selected_pdb_structural_template = {}
		if 'intact' in mutation_file:
			best_blast_alignment = pickle_load(osp.join(self.data_dir, 'interactome', "intact_blast_best_alignments_physical.pickle"))
			interfacial_residues = pickle_load(osp.join(self.data_dir, 'interactome', 'intact_interfacial_residues_for_selected_template.pickle'))
			selected_pdb_structural_template = pickle_load(osp.join(self.data_dir, 'interactome', 'intact_selected_pdb_structural_template.pickle'))
		else:
			best_blast_alignment = pickle_load(osp.join(self.data_dir, 'interactome', "hiunion_blast_best_alignments.pickle"))
			interfacial_residues = pickle_load(osp.join(self.data_dir, 'interactome', 'hiunion_interfacial_residues_for_selected_template.pickle'))
			selected_pdb_structural_template = pickle_load(osp.join(self.data_dir, 'interactome', 'hiunion_selected_pdb_structural_template.pickle'))

		# info to update
		interactors, missense_info_with_protein_partners = self.get_interactors(header_dict, missense_info, interfacial_residues) # key = uniprot protein, value = list of protein partners
		missense_info_position_within_blast_alignment = []

		total_num_mutations = 0
		num_removed = 0
		num_kept = 0
		missense_info_position_within_blast_alignment.append(missense_info_with_protein_partners[0])
		print('-----' + mutation_file + '-----')
		header_dict['protein_partners'] = len(header_dict)
		missense_info_position_within_blast_alignment[0].append('protein_partners')
		header_dict['on_interfacial_res'] = len(header_dict)
		missense_info_position_within_blast_alignment[0].append('on_interfacial_res') # on interfaicial residue
		header_dict['ir_mutation'] = len(header_dict)
		missense_info_position_within_blast_alignment[0].append('ir_mutation') # 0/1 (1 = on interfacial residue, 0 = not on interfacial residue)

		num_IR, num_non_IR = 0, 0
		for i in range(1, len(missense_info_with_protein_partners)):
			total_num_mutations += 1
			uniprot_protein = missense_info_with_protein_partners[i][header_dict['uniprot_protein']]
			mut_pos = int(missense_info_with_protein_partners[i][header_dict['prot_res_pos_in_uniprot_protein']]) # mutation position in uniprot numbering
			partners = interactors[uniprot_protein]
			# check that mutation is within blast alignment of for each interaction, remove interaction if is not
			# if no interactions left, then remove mutation
			confirmed_partners = [] # keep only interactions where mutation lies within blast alignment
			perturbations = []
			for p in partners:
				residues = []
				pdb, chain1, chain2 = "", "", ""
				# get residues and perturbations
				if (uniprot_protein, p) in interfacial_residues:
					residues = interfacial_residues[(uniprot_protein, p)][0]
					template = selected_pdb_structural_template[(uniprot_protein, p)]
					pdb, chain1, chain2 = template.split('_')
				elif (p, uniprot_protein) in interfacial_residues:
					residues = interfacial_residues[(p, uniprot_protein)][1]
					template = selected_pdb_structural_template[(p, uniprot_protein)]
					pdb, chain2, chain1 = template.split('_')
				else:
					print('Something went wrong! Cannot find interfacial residues for:', uniprot_protein, p)
					return
				# check that mutation is within blast alignment for template1
				# do not need to check for template2
				template1 = '_'.join([pdb, chain1])
				alignment = best_blast_alignment[(uniprot_protein, template1)]
				qstart, qend = int(alignment[4]), int(alignment[5])
				if mut_pos >= qstart and mut_pos <= qend:
					confirmed_partners.append(p)
					if mut_pos in residues: # see if mutation position coincides with a residue
						perturbations.append('1')
					else:
						perturbations.append('0')
				# else:
				# 	print('not within blast alignment:', uniprot_protein, template1, mut_pos, qstart, qend)
			if len(confirmed_partners) != len(perturbations):
				print('Something went wrong!')
			else:
				if len(confirmed_partners) == 0:
					# print('no interactions left for mutation:', uniprot_protein, partners)
					num_removed += 1
				else:
					num_kept += 1
					# add protein partners, perturbations, and edgotype to final_missense_info
					missense_info_with_protein_partners[i].append(','.join(confirmed_partners))
					missense_info_with_protein_partners[i].append(','.join(perturbations))
					if '1' in perturbations: 
						missense_info_with_protein_partners[i].append('1')
						num_IR += 1
					else:
						missense_info_with_protein_partners[i].append('0')
						num_non_IR += 1
					missense_info_position_within_blast_alignment.append(missense_info_with_protein_partners[i])

		print('Total number of mutations with protein partners:', total_num_mutations)
		print('Number of mutations removed:', num_removed)
		print('Number of mutations kept (with mutation positions within blast alignments):', num_kept)
		print('Number of IR mutations:', num_IR, round(num_IR/num_kept*100, 2))
		print('Number of NIR (non-IR) mutations:', num_non_IR, round(num_non_IR/num_kept*100, 2))

		write_missense_info_to_file(missense_info_position_within_blast_alignment, ir_mutations_file)

	def get_ir_mutations_all(self):
		for mutation_file in self.mutation_files_list:
			mutation_file_items = mutation_file.split('_')
			ir_mutations_file = ''
			if 'resistant' in mutation_file:
				ir_mutations_file = osp.join(self.data_dir, 'edgotypes', '_'.join([mutation_file_items[0], mutation_file_items[2], 'resistant_missense_mutations_ir.tsv']))
			elif 'pathogenic' in mutation_file:
				ir_mutations_file = osp.join(self.data_dir, 'edgotypes', '_'.join([mutation_file_items[0], mutation_file_items[2], 'missense_mutations_pathogenic_ir.tsv']))
			elif 'cancer_driving' in mutation_file:
				ir_mutations_file = osp.join(self.data_dir, 'edgotypes', '_'.join([mutation_file_items[0], 'cancer_driving', 'missense_mutations_ir.tsv']))
			else:
				ir_mutations_file = osp.join(self.data_dir, 'edgotypes', '_'.join([mutation_file_items[0], mutation_file_items[2], 'missense_mutations_ir.tsv']))

			self.find_ir_mutations(mutation_file, ir_mutations_file)


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')
	final_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'mutations_final')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	dbsnp_pathogenicities = ['', 'pathogenic']
	mutation_files_list = []

	for interactome in interactomes:
		for data in mutation_data:
			if data == 'dbsnp':
				for dbsnp_pathogenicity in dbsnp_pathogenicities:
					if dbsnp_pathogenicity == '':
						mutation_files_list.append('_'.join([interactome, 'mapped', data, 'missense', 'mutations', 'nonredundant.tsv']))
					else:
						mutation_files_list.append('_'.join([interactome, 'mapped', data, 'missense', 'mutations', dbsnp_pathogenicity, 'nonredundant.tsv']))
			else: # data == 'clinvar'
				mutation_files_list.append('_'.join([interactome, 'mapped', data, 'missense', 'mutations', 'nonredundant.tsv']))

	i = IRMutation(data_dir, final_data_dir, mutation_files_list)
	i.get_ir_mutations_all()

if __name__=='__main__':
	main()