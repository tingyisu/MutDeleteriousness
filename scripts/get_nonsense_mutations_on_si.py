'''
Maps nonsense mutations to the HI-union and IntAct structural interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import os.path as osp
import os
from simple_tools import pickle_load, get_missense_info, write_missense_info_to_file

# removes mutations that do not have protein partners or that don't have any leftoever interactions in which the mutation position lies within the structural template
# finds where mutation lies on protein (first-half or last-half)

class NonsenseMutation:
	def __init__(self, data_dir, final_data_dir, mutation_files_list):
		self.data_dir = data_dir
		self.swissprot_seq_dict = pickle_load(osp.join(self.data_dir, 'swissprot_seq_dict.pickle'))
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


	# find nonsense mutations that fall within BLAST alignment
	def find_mutations_on_si(self, mutation_file, nonsense_mutations_file):
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
		num_in_first_half, num_in_last_half = 0, 0
		missense_info_position_within_blast_alignment.append(missense_info_with_protein_partners[0])
		print('-----' + mutation_file + '-----')
		header_dict['protein_partners'] = len(header_dict)
		missense_info_position_within_blast_alignment[0].append('protein_partners')

		for i in range(1, len(missense_info_with_protein_partners)):
			total_num_mutations += 1
			uniprot_protein = missense_info_with_protein_partners[i][header_dict['uniprot_protein']]
			mut_pos = int(missense_info_with_protein_partners[i][header_dict['prot_res_pos_in_uniprot_protein']]) # mutation position in uniprot numbering
			partners = interactors[uniprot_protein]
			# check that mutation is within blast alignment of for each interaction, remove interaction if is not
			# if no interactions left, then remove mutation
			confirmed_partners = [] # keep only interactions where mutation lies within blast alignment
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
				# else:
				# 	print('not within blast alignment:', uniprot_protein, template1, mut_pos, qstart, qend)


			# make sure that mutation lies in at least one blast alignment
			# and that it causes a truncation of >= 50% of the protein
			uniprot_protein = missense_info_with_protein_partners[i][header_dict['uniprot_protein']]
			prot_res_pos_in_uniprot_protein = int(missense_info_with_protein_partners[i][header_dict['prot_res_pos_in_uniprot_protein']])
			len_of_uniprot_protein = len(self.swissprot_seq_dict[uniprot_protein])
			# print(prot_res_pos_in_uniprot_protein, len_of_uniprot_protein)


			if len(confirmed_partners) == 0:
				# print('no interactions left for mutation:', uniprot_protein, partners)
				num_removed += 1
			else:
				if prot_res_pos_in_uniprot_protein/len_of_uniprot_protein <= 0.5:
					num_in_first_half += 1
				else:
					num_in_last_half += 1


				num_kept += 1
				# add protein partners, perturbations, and edgotype to final_missense_info
				missense_info_with_protein_partners[i].append(','.join(confirmed_partners))
				missense_info_position_within_blast_alignment.append(missense_info_with_protein_partners[i])

		print('Total number of mutations with protein partners:', total_num_mutations)
		print('Number of mutations removed:', num_removed)
		# print('Number of mutations kept (with mutation positions within blast alignments & truncation of >= 50%):', num_kept)
		print('Number of mutations kept (with mutation positions within blast alignments):', num_kept)
		print('Number of mutations in first & last halves:', num_in_first_half, num_in_last_half)


		write_missense_info_to_file(missense_info_position_within_blast_alignment, nonsense_mutations_file)

	def get_nonsense_mutations_all(self):
		for mutation_file in self.mutation_files_list:
			# print('-----' + mutation_file + '-----')
			mutation_file_items = mutation_file.split('_')
			nonsense_mutations_file = osp.join(self.data_dir, 'nonsense_on_si', '_'.join([mutation_file_items[0], mutation_file_items[2], 'nonsense_mutations.tsv']))
			# print(nonsense_mutations_file)
			self.find_mutations_on_si(mutation_file, nonsense_mutations_file)


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')
	final_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'mutations_final')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	mutation_files_list = []

	for interactome in interactomes:
		for data in mutation_data:
			mutation_files_list.append('_'.join([interactome, 'mapped', data, 'nonsense', 'mutations', 'nonredundant.tsv']))

	i = NonsenseMutation(data_dir, final_data_dir, mutation_files_list)
	i.get_nonsense_mutations_all()


if __name__=='__main__':
	main()