'''
Processes binary interactions in the HI-union dataset of human protein-protein interactions (http://www.interactome-atlas.org/download)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from Bio import SeqIO
import pickle
import pandas as pd
from simple_tools import get_list_from_file, pickle_dump, pickle_load

class Data: 
	def __init__(self, orig_data_dir, processed_data_dir):
		# files to process
		self.ref_ppi_file = osp.join(orig_data_dir, 'HI-union.tsv')
		self.pdb_seqres_file = osp.join(orig_data_dir, 'pdb_seqres.txt')
		# attributes to load/find
		self.swissprot_seq_dict = pickle_load(osp.join(processed_data_dir, 'swissprot_seq_dict.pickle'))
		self.ensembl_uniprot = pickle_load(osp.join(processed_data_dir, 'ensembl_gene_uniprot_dict.pickle')) # dict with key = ensembl gene id, value = list of swissprot (uniprot) ids
		self.pairs = [] # list of tuples of interacting Ensembl IDs and their corresponding Uniprot IDs (Ensembl1, Ensembl2, Uniprot1, Uniprot2)
		self.unique_uniprot = [] # list of unique uniprot proteins
		self.unique_uniprot_pairs = [] # list of unique uniprot interacting pairs
		# files to write to
		interactome_data_dir = osp.join(processed_data_dir, 'interactome')
		self.pairs_pickle_file = osp.join(interactome_data_dir, 'HI-union_interacting_pairs.pickle')
		self.unique_uniprot_pairs_pickle_file = osp.join(interactome_data_dir, 'HI-union_unique_uniprot_pairs.pickle')
		self.uniprot_fasta_blast_file = osp.join(interactome_data_dir, 'hiunion_uniprot_sequences_blast.fasta')
		self.pdb_seqres_blast_file = osp.join(interactome_data_dir, 'pdb_seqres_blast.fasta')

	# updates self.pairs, self.unique_uniprot, and self.unique_uniprot_pairs
	def get_data(self):
		print('------Getting a list of interacting Ensembl/Uniprot IDs and unique UniProt proteins/pairs-----')
		num_self = 0
		with open(self.ref_ppi_file, 'r') as fname:
			for line in fname:
				ensembl1, ensembl2 = line.strip('\n').split('\t')

				# need both interaction Ensembl genes to have UniProtKB/SWISS-PROT reviewed proteins
				if ensembl1 in self.ensembl_uniprot and ensembl2 in self.ensembl_uniprot: # some ensembl ids do not have swissprot ids

					# add to self.pairs
					uniprot1, uniprot2 = self.ensembl_uniprot[ensembl1], self.ensembl_uniprot[ensembl2] 

					if uniprot1 != uniprot2:
						# may have multiple Ensembl pairs with the same UniProt pair...or the reverse...just pick the first one encountered
						if (uniprot1, uniprot2) not in self.unique_uniprot_pairs and (uniprot2, uniprot1) not in self.unique_uniprot_pairs:
							self.unique_uniprot_pairs.append((uniprot1, uniprot2))
							# add to self.unique_uniprot
							self.unique_uniprot.append(uniprot1)
							self.unique_uniprot.append(uniprot2)
							self.pairs.append((ensembl1, ensembl2, uniprot1, uniprot2))
						# else:
						# 	print('repeated:', (uniprot1, uniprot2))
					else:
						num_self += 1

		self.unique_uniprot = list(set(self.unique_uniprot))


		print('Number of self-interactions (removed):', num_self)
		print('Number of Ensembl IDs with corresponding Swissprot/Uniprot IDs:', len(self.ensembl_uniprot))
		print('Number of PPIs in HI-union: 64006')
		print('Number of PPIs found:',  len(self.pairs))
		print('Number of unique Uniprot pairs:', len(self.unique_uniprot_pairs))
		print('Number of unique Uniprot proteins:', len(self.unique_uniprot))

		# pickles self.uniprot_pairs to be used in build_structural_interactome.py
		pickle_dump(self.pairs, self.pairs_pickle_file)

		# save to pickle file
		pickle_dump(self.unique_uniprot_pairs, self.unique_uniprot_pairs_pickle_file)


	# get file for blast
	def get_uniprot_fasta_blast(self):
		print('-----Writing Uniprot Seq to fasta BLAST format-----')
		with open(self.uniprot_fasta_blast_file, 'w') as f:
			for p in self.unique_uniprot:
				f.write('>' + str(p) + '\n')
				f.write(str(self.swissprot_seq_dict[p]) + '\n')

	# converts pdb_seqres.txt file (downloaded from PDB) to BLAST usable form
	def to_blast(self):
		print('-----Converting pdb_seqres.txt to BLAST format-----')
		with open(self.pdb_seqres_file, 'r') as fasta, open(self.pdb_seqres_blast_file, 'w') as to_blast:
			i = 0
			name = ""
			for line in fasta.readlines():
				i += 1
				if i % 2 == 1:
					name = line.split(':')[0]
					name = name.split(' ')[0]
					to_blast.write(name + "\n")
				else: to_blast.write(line)

	def process_hiunion_data(self):
		self.get_data()
		self.get_uniprot_fasta_blast()
		self.to_blast()

def main():
	script_dir = osp.dirname(__file__)
	orig_data_dir = osp.join(script_dir, '..', 'data', 'original')
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')
	
	d = Data(orig_data_dir, processed_data_dir)
	d.process_hiunion_data()

if __name__=='__main__':
	main()