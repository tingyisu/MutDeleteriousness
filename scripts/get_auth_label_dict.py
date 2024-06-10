'''
Creates a dictionary of PDB auth_seq_id (author-specified numbering) to label_seq_id (PDB assigned sequential numbering) mappings for residues
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import argparse
import os.path as osp
from Bio.PDB import PDBList
import os
from simple_tools import download_pdb_structure, pickle_dump

class AuthLabel:
	def __init__(self, name, pdb_download_path, output_directory):
		self.pdb_download_path = pdb_download_path
		self.name = name
		if self.name == 'intact':
			self.interactions_file = osp.join(output_directory, self.name + "_interfacial_residues_physical.tsv")
		else:
			self.interactions_file = osp.join(output_directory, self.name + "_interfacial_residues.tsv")
		self.auth_label_dict_pickle_file = osp.join(output_directory, self.name + "_auth_label_dict.pickle")
		self.uniprot_pdb_chains = {} # key = uniprot protein, value = list of pdb chains (confirmed in interactions_file)
		self.uniprot_pdb_chains_pickle_file = osp.join(output_directory, self.name + '_uniprot_pdb_chains.pickle')
		self.pdb_chains = {} # key = pdb structure, value = list of chains (confirmed in interactions_file)
		self.auth_label_dict = {} # dict of dict of dict of key = auth_seq_id, value = label_seq_id, (PDB - pdb_seq_res.txt and biopython both use auth as convention)

	def get_uniprot_pdb_chains(self):
		with open(self.interactions_file, 'r') as f:
			next(f)
			for line in f:
				items = line.strip().split('\t')
				uniprot1, uniprot2 = "", ""
				pdb_chain1, pdb_chain2 = "", ""
				if self.name == 'hiunion':
					uniprot1, uniprot2 = items[2], items[3]
					pdb_chain1, pdb_chain2 = items[4] + '_' + items[5], items[4] + '_' + items[6]
				else:
					uniprot1, uniprot2 = items[4], items[5]
					pdb_chain1, pdb_chain2 = items[6] + '_' + items[7], items[6] + '_' + items[8]
				self.uniprot_pdb_chains.setdefault(uniprot1, []).append(pdb_chain1)
				self.uniprot_pdb_chains.setdefault(uniprot2, []).append(pdb_chain2)

		# get unique pdb+chains only
		for uniprot in self.uniprot_pdb_chains:
			unique_chains = list(set(self.uniprot_pdb_chains[uniprot]))
			self.uniprot_pdb_chains[uniprot] = unique_chains
		pickle_dump(self.uniprot_pdb_chains, self.uniprot_pdb_chains_pickle_file)

	def get_pdb_chains(self):

		for uniprot in self.uniprot_pdb_chains:
			for pdb_chain in self.uniprot_pdb_chains[uniprot]:
				pdb, chain = pdb_chain.split('_')
				self.pdb_chains.setdefault(pdb, []).append(chain)

		# get unique chains only
		for pdb in self.pdb_chains:
			unique_chains = list(set(self.pdb_chains[pdb]))
			self.pdb_chains[pdb] = unique_chains
		# print(self.pdb_chains['4awl'])

	# gets self.auth_label_dict
	# since PDB & Biopython both use auth (auth_seq_id) instead of the actual position on sequence, need to get the corresponding label_seq_id
	# necessary for mapping residues to uniprot proteins in BLASTP alignments
	'''
	format of self.auth_label_dict 
	{
	'4awl' :
			{
			'A': {14: 24, ..., 89: 99}
			}
	}
	'''
	def get_auth_label_dict(self):
		print('-----Getting self.auth_label_dict-----')
		print('There are', len(self.pdb_chains), 'PDB structures')
		i = 0
		for pdb in self.pdb_chains:
			has_column_names = False
			if not osp.exists(osp.join(self.pdb_download_path, pdb + '.cif')):
				download_pdb_structure(pdb, self.pdb_download_path)
			i += 1
			fname = osp.join(self.pdb_download_path, pdb + '.cif')
			self.auth_label_dict[pdb] = {}
			column_names = []
			with open(fname, 'r') as f:
				for line in f:
					# get column_names (e.g. ['_atom_side.PDB_group', '_atom_site.id', ..., '_atom_site.pdbx_PDB_model_num'])
					if line.startswith('_atom_site.'):
						column_names.append(line.strip())
						if not has_column_names:
							has_column_names = True

					elif line.startswith('ATOM') and not line.startswith('ATOMS') and has_column_names:
						items = line.strip().split()
						if len(column_names) != len(items):
							print('ERROR! column lengths do not match! ', pdb, items)

						# get indices
						chain_index = column_names.index('_atom_site.auth_asym_id')
						auth_res_index = column_names.index('_atom_site.auth_seq_id')
						auth_res_identifier_index = column_names.index('_atom_site.pdbx_PDB_ins_code') # some residues have the same auth_seq_id but a different 'identifier'/PDB ins code (e.g. 1kmc_B has 175, 175_A, 175_B, 175_C)
						label_res_index = column_names.index('_atom_site.label_seq_id')

						# get chain id and residue ids (in both auth and label format)
						chain, auth_res, auth_res_identifier, label_res = items[chain_index], items[auth_res_index], items[auth_res_identifier_index], items[label_res_index]

						# get complete auth_seq_id along with the identifier
						auth_res_complete = auth_res
						if auth_res_identifier != '?':
							auth_res_complete = auth_res + auth_res_identifier

						# create dictionary for chain in self.auth_label_dict
						if chain not in self.auth_label_dict[pdb]:
							self.auth_label_dict[pdb][chain] = {}

						# add complete auth_seq_id with the identifier
						if auth_res_complete not in self.auth_label_dict[pdb][chain]:
							self.auth_label_dict[pdb][chain][auth_res_complete] = label_res
						else: 
							prev_label_res = self.auth_label_dict[pdb][chain][auth_res_complete]
							if prev_label_res != label_res:
								print('Have the same auth_seq_id for different label_seq_ids!', pdb, chain, auth_res_complete, prev_label_res, label_res)

					elif line.startswith('_pdbx_poly_seq_scheme.') or line.startswith('_atom_site_anisotrop.'): # end of atoms (and hetatoms) lines
						break
					else:
						continue
			print('Done looking at', i, '/', len(self.pdb_chains), 'structures')

		pickle_dump(self.auth_label_dict, self.auth_label_dict_pickle_file)

	def get_chains_auth_label_dict(self):
		self.get_uniprot_pdb_chains()
		self.get_pdb_chains()
		self.get_auth_label_dict()

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--pdb_download_path') # e.g. /home/username/scratch/pdb_cif
	parser.add_argument('-o', '--output_directory')  # e.g. /home/username/scratch/interfacial_residues/files
	args = parser.parse_args()

	hiunion = AuthLabel('hiunion', args.pdb_download_path, args.output_directory)
	hiunion.get_chains_auth_label_dict()

	intact = AuthLabel('intact', args.pdb_download_path, args.output_directory)
	intact.get_chains_auth_label_dict()


if __name__=='__main__':
	main()