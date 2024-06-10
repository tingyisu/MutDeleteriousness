'''
Combine BLASTP alignments from HI-Union and IntAct interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import argparse
import pickle
from simple_tools import pickle_load, pickle_dump, check_exist, download_pdb_structure

class Convert:
	def __init__(self, data_dir, pdb_download_dir):
		self.data_dir = data_dir
		self.blast_alignments = pickle_load(osp.join(self.data_dir, 'all_blast_best_alignments_reduced_info.pickle')) # key = (protein, pdb_chain), value = [pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, protein alignment, template_chain alignment]
		self.modeller_templates = pickle_load(osp.join(self.data_dir, 'all_modeller_templates.pickle'))
		self.pdb_download_dir = pdb_download_dir
		self.pdb_structures = pickle_load(osp.join(self.data_dir, 'pdb_structures_to_download.pickle'))
		self.auth_label_dict = {}
		self.auth_label_dict_pickle_file = osp.join(self.data_dir, 'auth_label_dict_modeller.pickle') # dict of dict of dict to convert auth_seq_id to label_seq_id
		self.label_auth_dict = {} # dict of dict of dict to convert label_seq_id to auth_seq_id
		self.label_auth_dict_pickle_file = osp.join(self.data_dir, 'label_auth_dict_modeller.pickle')
		# self.auth_label_dict = pickle_load(osp.join(self.data_dir, self.name + '_auth_label_dict_modeller.pickle'))
		# self.label_auth_dict = pickle_load(osp.join(self.data_dir, self.name + '_label_auth_dict_modeller.pickle'))
		self.converted_blast_alignments = {} # converted so that the BLAST alignment has gaps for residues that don't appear in the mmCIF (PDB) file (on the subject side)
		self.converted_blast_alignments_pickle_file = osp.join(self.data_dir, 'converted_blast_alignments.pickle')
		# need for renumbering residues back to original uniprot alignment indices (on the query side)
		self.uniprot_alignment_dict = {} # key = (uniprot, pdb_chain), value = dict of modeller residue numbering (starting from one), and original uniprot residue numbering key-value pairs
		self.uniprot_alignment_dict_pickle_file = osp.join(self.data_dir, 'uniprot_alignment_dict.pickle')


	def get_auth_label_dict(self):
		print('-----Getting self.auth_label_dict-----')
		print('There are', len(self.pdb_structures), 'PDB structures')
		i = 0
		for pdb in self.pdb_structures:
			has_column_names = False
			if not osp.exists(osp.join(self.pdb_download_dir, pdb + '.cif')):
				download_pdb_structure(pdb, self.pdb_download_dir)
			i += 1
			fname = osp.join(self.pdb_download_dir, pdb + '.cif')
			self.auth_label_dict[pdb] = {}
			column_names = []
			with open(fname, 'r') as f:
				for line in f:
					# get column_names (e.g. ['_atom_side.PDB_group', '_atom_site.id', ..., '_atom_site.pdbx_PDB_model_num'])
					if line.startswith('_atom_site.'):
						column_names.append(line.strip())
						if not has_column_names:
							has_column_names = True

					elif (line.startswith('ATOM') or line.startswith('HETATM')) and not line.startswith('ATOMS') and has_column_names:
						items = line.split() # split by whitespaces, not '\t'! splitting by whitepaces removes \n and splits by \t as well!
						if len(column_names) != len(items):
							print('ERROR! column lengths do not match!', pdb, items)

						# get indices
						residue_index = column_names.index('_atom_site.label_comp_id') # residue (in 3 letter code), should be the same as the auth_comp_id
						residue_index2 = column_names.index('_atom_site.auth_comp_id') # residue (in 3 letter code)
						chain_index = column_names.index('_atom_site.auth_asym_id')
						auth_res_index = column_names.index('_atom_site.auth_seq_id')
						auth_res_identifier_index = column_names.index('_atom_site.pdbx_PDB_ins_code') # some residues have the same auth_seq_id but a different 'identifier'/PDB ins code (e.g. 1kmc_B has 175, 175_A, 175_B, 175_C)
						label_res_index = column_names.index('_atom_site.label_seq_id')


						# get chain id and residue ids (in both auth and label format)
						residue, residue2, chain, auth_res, auth_res_identifier, label_res = items[residue_index], items[residue_index2], items[chain_index], items[auth_res_index], items[auth_res_identifier_index], items[label_res_index]

						if residue != residue2:
							print('ERROR! Different residues in label_comp_id and auth_comp_id', fname, residue, residue2)
						else:
							if (line.startswith('HETATM') and residue == 'MSE') or line.startswith('ATOM'): # get HETATM if is MSE https://salilab.org/archives/modeller_usage/2009/msg00185.html
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
			print('Done looking at', i, '/', len(self.pdb_structures), 'structures')
			# print('Done checking', i, '/', len(self.pdb_chains), 'structures')
		pickle_dump(self.auth_label_dict, self.auth_label_dict_pickle_file)


	def get_label_auth_dict(self):
		# reverse auth_label_dict to get label_auth_dict
		# if not check_exist(self.label_auth_dict_pickle_file):
		for pdb_struct in self.auth_label_dict:
			self.label_auth_dict[pdb_struct] = {}
			for chain in self.auth_label_dict[pdb_struct]:
				self.label_auth_dict[pdb_struct][chain] = {}
				for auth_id in self.auth_label_dict[pdb_struct][chain]:
					label_id = self.auth_label_dict[pdb_struct][chain][auth_id]
					self.label_auth_dict[pdb_struct][chain][label_id] = auth_id
		pickle_dump(self.label_auth_dict, self.label_auth_dict_pickle_file)


	def convert_blast_alignments(self):
		# convert sstart and send from label_seq_id format to auth_seq_id format
		for modeller_template in self.modeller_templates:
			uniprot, pdb, chain = modeller_template.split('_')
			pdb_chain = '_'.join([pdb, chain])
			qstart, qend, label_sstart, label_send, qalignment, salignment = self.blast_alignments[(uniprot, pdb_chain)]
			# get only label_seq_ids that exist in the mmCIF file
			label_residues_that_exist = [str(i) for i in range(int(label_sstart), int(label_send) + 1) if str(i) in self.label_auth_dict[pdb][chain]] # label_seq_ids that exist in the mmCIF file
			converted_salignment = ''
			# change residues in salignment that do not exist in the mmCIF file to '-'
			i = int(label_sstart)
			for residue in salignment:
				if residue != '-':
					if str(i) in label_residues_that_exist:
						converted_salignment += residue
					else:
						converted_salignment += '-'
					i += 1
				else:
					converted_salignment += '-'
			if i != (int(label_send) + 1):
				print('ERROR! Something went wrong with the label iterations!!', i, label_send)
			auth_sstart = self.label_auth_dict[pdb][chain][label_residues_that_exist[0]]
			auth_send = self.label_auth_dict[pdb][chain][label_residues_that_exist[len(label_residues_that_exist)-1]]
			# get uniprot alignment dict
			j = int(qstart)
			curr_index = 1
			self.uniprot_alignment_dict['_'.join([uniprot, pdb_chain])] = {}
			for residue in qalignment:
				if residue != '-':
					self.uniprot_alignment_dict['_'.join([uniprot, pdb_chain])][curr_index] = j
					j += 1
					curr_index += 1
				else:
					continue
			if j != (int(qend) + 1):
				print('ERROR! Something went wrong with the query (uniprot) iterations!!', j, qend)

			self.converted_blast_alignments[(uniprot, pdb_chain)] = [qstart, qend, str(auth_sstart), str(auth_send), qalignment, converted_salignment]


		pickle_dump(self.converted_blast_alignments, self.converted_blast_alignments_pickle_file)
		pickle_dump(self.uniprot_alignment_dict, self.uniprot_alignment_dict_pickle_file)

	def convert_blast_alignment_to_pdb_numbering(self):
		self.get_auth_label_dict()
		self.get_label_auth_dict()
		self.convert_blast_alignments()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data_dir')  #/home/username/projects/def-*/username/modeller/files
	parser.add_argument('-p', '--pdb_download_dir') # /home/username/scratch/pdb_cif
	args = parser.parse_args()

	print('Converting BLAST alignments...')
	c = Convert(args.data_dir, args.pdb_download_dir)
	c.convert_blast_alignment_to_pdb_numbering()


if __name__=='__main__':
	main()