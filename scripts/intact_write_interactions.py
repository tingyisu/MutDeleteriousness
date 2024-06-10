'''
Writes interfacial residues for IntAct PDB chain-pairs to file
	# 1. get memoized interfacial residues from file in format (chain_pair1, chain_pair2) = [res1_list, res2_list]
	# 2. get interactions dict
	# 3. write interactions dict to file
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import argparse
import os.path as osp
from simple_tools import pickle_load

class InteractionsFile:
	def __init__(self, pos_hom_file, gene_dict_file, intact_id_dict_file, memo_residues_file, interactions_file):
		self.pos_hom = pickle_load(pos_hom_file)
		self.gene_dict = pickle_load(gene_dict_file) # key = uniprot protein id, value = gene name
		self.intact_id_dict = pickle_load(intact_id_dict_file) # key = uniprot protein id, value = intact id
		self.memo_residues_file = memo_residues_file
		self.interactions_file = interactions_file
		self.memo_residues = {}

	def get_memo_from_file(self):
		print('-----Reading in residues from memo file-----')
		with open(self.memo_residues_file, 'r') as fname:
			for line in fname:
				if line.strip(): # if line is not empty
					items = line.strip().split('\t')
					struct_name, chain1_name, chain2_name = items[0], items[1], items[2]
					chain1, chain2 = struct_name + '_' + chain1_name, struct_name + '_' + chain2_name
					if len(items) == 3: # no interfacial residues
						self.memo_residues[(chain1, chain2)] = [[],[]]
					else:
						residue1_str, residue2_str = items[3], items[4]
						self.memo_residues[(chain1, chain2)] = [residue1_str.split(','), residue2_str.split(',')]
		print('Number of memoized chain pairs:', len(self.memo_residues))

	# get interactions and write them to file
	def get_interactions(self, inc_self, self_only):
		print('-----Getting and writing interactions to file-----')
		uniprot_pairs = []
		with open(self.interactions_file, 'w') as f:
			f.write('\t'.join(['Gene1_Name', 'Gene2_Name', 'Protein1_IntAct_ID', 'Protein2_IntAct_ID', 'Protein1', 'Protein2', 'PDB_Structure', 'Protein1_Chain', 'Protein2_Chain', 'Protein1_Interfacial_Residues', 'Protein2_Interfacial_Residues']) + '\n')
			for pair in self.pos_hom:
				p1, p2 = pair
				if not inc_self and p1 == p2:
					continue

				# if want self interactions only
				if self_only and p1 != p2:
					continue

				chain_pairs_list = self.pos_hom[(p1, p2)]
				for (chain1, chain2) in chain_pairs_list:
					if (chain1, chain2) in self.memo_residues:
						res1_list, res2_list = self.memo_residues[(chain1, chain2)]
						if res1_list != [] and res2_list != []:
							if pair not in uniprot_pairs:
									uniprot_pairs.append(pair)
							pdb_struct, chain1_name, chain2_name = chain1.split('_')[0], chain1.split('_')[1], chain2.split('_')[1]
							f.write('\t'.join([self.gene_dict[p1], self.gene_dict[p2], self.intact_id_dict[p1], self.intact_id_dict[p2], p1, p2, pdb_struct, chain1_name, chain2_name, ','.join(res1_list), ','.join(res2_list)]) + '\n')
					elif (chain2, chain1) in self.memo_residues:
						res2_list, res1_list = self.memo_residues[(chain2, chain1)]
						if res2_list != [] and res1_list != []:
							if pair not in uniprot_pairs:
									uniprot_pairs.append(pair)
							pdb_struct, chain2_name, chain1_name = chain2.split('_')[0], chain2.split('_')[1], chain1.split('_')[1]
							f.write('\t'.join([self.gene_dict[p1], self.gene_dict[p2], self.intact_id_dict[p1], self.intact_id_dict[p2], p1, p2, pdb_struct, chain1_name, chain2_name, ','.join(res1_list), ','.join(res2_list)]) + '\n')
					else:
						print('Missing!!', (chain1, chain2))
						
		print('Number of PPIs:', len(uniprot_pairs))
		for (p1, p2) in uniprot_pairs:
			if (p2, p1) in uniprot_pairs:
				print('repeated!', (p1, p2), (p2, p1))
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-o', '--output_directory') # /home/username/scratch/interfacial_residues/files
	args = parser.parse_args()
	pos_hom_file = osp.join(args.output_directory,"intact_pos_hom_physical.pickle")
	gene_dict_file = osp.join(args.output_directory,"intact_gene_dict_physical.pickle")
	intact_id_dict_file = osp.join(args.output_directory,"intact_id_dict_physical.pickle")
	combined_memo_residues_file = osp.join(args.output_directory, "memo_residues_combined.tsv") # hiunion + intact (since have some overlapping)
	interactions_file = osp.join(args.output_directory, "intact_interfacial_residues_physical.tsv")

	intact = InteractionsFile(pos_hom_file, gene_dict_file, intact_id_dict_file, combined_memo_residues_file, interactions_file)
	intact.get_memo_from_file()
	intact.get_interactions(False, False)


if __name__=='__main__':
	main()
import pickle
