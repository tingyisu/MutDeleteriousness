'''
Writes interfacial residues for HI-union PDB chain-pairs to file
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
	def __init__(self, pairs_pickle_file, pos_hom_file, memo_residues_file, interactions_file):
		self.pairs = pickle_load(pairs_pickle_file) # list of tuples of interacting Ensembl IDs and their corresponding Uniprot IDs (Ensembl1, Ensembl2, Uniprot1, Uniprot2)
		# self.unique_uniprot_pairs = [] # list of unique uniprot interacting pairs
		# self.uniprot_to_ensembl = {} # key = (Uniprot1, Uniprot2), value = (Ensembl1, Ensembl2)
		self.pos_hom = pickle_load(pos_hom_file)
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
			f.write('\t'.join(['Ensembl_Gene_ID1', 'Ensembl_Gene_ID2', 'Protein1', 'Protein2', 'PDB_Structure', 'Protein1_Chain', 'Protein2_Chain', 'Protein1_Interfacial_Residues', 'Protein2_Interfacial_Residues']) + '\n')
			for pair in self.pairs:
				e1, e2, p1, p2 = pair
				if not inc_self and p1 == p2: # probably don't need this anymore, since already removed self-interactions in process_hiunion_data.py
					continue

				# if want self interactions only
				if self_only and p1 != p2:
					continue

				if (p1, p2) in self.pos_hom: # may not be if have no possible interacting chains in any PDB structure
					chain_pairs_list = self.pos_hom[(p1, p2)]
					# print((p1, p2))
					for (chain1, chain2) in chain_pairs_list:
						if (chain1, chain2) in self.memo_residues:
							res1_list, res2_list = self.memo_residues[(chain1, chain2)]
							if res1_list != [] and res2_list != []:
								if (p1, p2) not in uniprot_pairs:
									uniprot_pairs.append((p1, p2))
								pdb_struct, chain1_name, chain2_name = chain1.split('_')[0], chain1.split('_')[1], chain2.split('_')[1]
								f.write('\t'.join([e1, e2, p1, p2, pdb_struct, chain1_name, chain2_name, ','.join(res1_list), ','.join(res2_list)]) + '\n')
						elif (chain2, chain1) in self.memo_residues:
							res2_list, res1_list = self.memo_residues[(chain2, chain1)]
							if res2_list != [] and res1_list != []:
								if (p1, p2) not in uniprot_pairs:
									uniprot_pairs.append((p1, p2))
								pdb_struct, chain2_name, chain1_name = chain2.split('_')[0], chain2.split('_')[1], chain1.split('_')[1]
								f.write('\t'.join([e1, e2, p1, p2, pdb_struct, chain1_name, chain2_name, ','.join(res1_list), ','.join(res2_list)]) + '\n')
						else:
							print('Missing!!', (chain1, chain2))

		unique_proteins = []
		print('Number of PPIs:', len(uniprot_pairs))
		for (p1, p2) in uniprot_pairs:
			if (p2, p1) in uniprot_pairs:
				print('repeated!', (p1, p2), (p2, p1))
			unique_proteins.extend([p1, p2])

		print('Number of unique proteins:', len(set(unique_proteins)))


	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-o', '--output_directory') # /home/username/scratch/interfacial_residues/files
	args = parser.parse_args()

	pairs_pickle_file = osp.join(args.output_directory, 'HI-union_interacting_pairs.pickle')
	pos_hom_file = osp.join(args.output_directory,"hiunion_pos_hom.pickle")
	combined_memo_residues_file = osp.join(args.output_directory, "hiunion_memo_residues_combined.tsv")
	interactions_file = osp.join(args.output_directory, "hiunion_interfacial_residues.tsv")  


	hiunion = InteractionsFile(pairs_pickle_file, pos_hom_file, combined_memo_residues_file, interactions_file)
	hiunion.get_memo_from_file()
	hiunion.get_interactions(False, False)
	


if __name__=='__main__':
	main()
