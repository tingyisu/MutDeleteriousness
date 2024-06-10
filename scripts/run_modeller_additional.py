'''
Run MODELLER on to build additional models needed for FoldX BuildModel (calculating mutation induced change on protein stability)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse
import os.path as osp 
from modeller import *
from modeller.automodel import *
from simple_tools import pickle_load, download_pdb_structure, pickle_dump, check_exist

class AdditionalModels:
	def __init__(self, script_dir, scratch_dir, pdb_download_dir):
		self.foldx_buildmodel_dir = osp.join(scratch_dir, 'foldx_buildmodel_all')
		self.data_dir = osp.join(script_dir, 'files')
		self.alignments_dir = osp.join(script_dir, 'alignments')
		self.results_dir = osp.join(script_dir, 'results')
		self.pdb_download_dir = pdb_download_dir
		self.blast_alignments = pickle_load(osp.join(self.data_dir, 'all_blast_best_alignments_reduced_info.pickle'))
		self.auth_label_dict = {}
		self.auth_label_dict_pickle_file = osp.join(self.data_dir, 'auth_label_dict_modeller_additional.pickle') # dict of dict of dict to convert auth_seq_id to label_seq_id
		self.label_auth_dict = {} # dict of dict of dict to convert label_seq_id to auth_seq_id
		self.label_auth_dict_pickle_file = osp.join(self.data_dir, 'label_auth_dict_modeller_additional.pickle')
		self.additional_modeller_templates = pickle_load(osp.join(self.foldx_buildmodel_dir, 'additional_modeller_templates.pickle'))
		self.pdb_structures = list(set([template.split('_')[1] for template in self.additional_modeller_templates]))
		print('Number of unique PDB structures:', len(self.pdb_structures))
		self.converted_blast_alignments = {}
		self.converted_blast_alignments_pickle_file = osp.join(self.data_dir, 'converted_blast_alignments_additional.pickle')
		self.uniprot_alignment_dict = {}
		self.uniprot_alignment_dict_pickle_file = osp.join(self.data_dir, 'uniprot_alignment_dict_additional.pickle')

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
						residue_index = column_names.index('_atom_site.label_comp_id')
						chain_index = column_names.index('_atom_site.auth_asym_id')
						auth_res_index = column_names.index('_atom_site.auth_seq_id')
						auth_res_identifier_index = column_names.index('_atom_site.pdbx_PDB_ins_code') # some residues have the same auth_seq_id but a different 'identifier'/PDB ins code (e.g. 1kmc_B has 175, 175_A, 175_B, 175_C)
						label_res_index = column_names.index('_atom_site.label_seq_id')

						# get chain id and residue ids (in both auth and label format)
						residue, chain, auth_res, auth_res_identifier, label_res = items[residue_index], items[chain_index], items[auth_res_index], items[auth_res_identifier_index], items[label_res_index]

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
		for modeller_template in self.additional_modeller_templates:
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

	def create_uniprot_ali_files(self): # put uniprot sequence into .ali format specified by MODELLER
		for (uniprot, pdb_chain) in self.converted_blast_alignments:
			pdb, chain = pdb_chain.split('_')
			fname = '_'.join([uniprot, pdb_chain])
			qstart, qend, auth_sstart, auth_send, qalignment, converted_salignment = self.converted_blast_alignments[(uniprot, pdb_chain)]
			print('Creating .ali file for', uniprot, pdb_chain + '...')
			fname = osp.join(self.alignments_dir, fname + '.ali')
			cif_fname = osp.join(self.pdb_download_dir, pdb + '.cif')
			# if not check_exist(fname):
			with open(fname, 'w') as f:
				f.write('>P1;' + pdb_chain + '\n')
				f.write('structure:' + cif_fname + ':' + auth_sstart + ':' + chain + ':' + auth_send + ':' + chain + '::::' + '\n')
				f.write(converted_salignment + '*' + '\n' + '\n')
				f.write('>P1;' + uniprot + '\n')
				f.write('sequence:' + uniprot + ':' + qstart + '::' + qend + ':::::' '\n')
				f.write(qalignment + '*' + '\n')

	# run this for models with initial alignment discrepancies
	def build_additional_models(self):
		for additional_modeller_template in self.additional_modeller_templates:
			uniprot, pdb, chain = additional_modeller_template.split('_')
			pdb_chain = '_'.join([pdb, chain])
			print('Finding model for:', uniprot, pdb_chain + '...')
			# run if haven't found model already
			model_name = osp.join(self.results_dir, additional_modeller_template + '.B99990001.pdb')
			if not check_exist(model_name):
				try:
					fname = '_'.join([uniprot, pdb_chain])
					log.minimal()
					env = Environ()
					env.io.hetatm = False # (NOTE: MODELLER seems to take HETATM MSE regardless of whether this is set to True/False)! It converts MSE to MET! https://salilab.org/archives/modeller_usage/2009/msg00185.html
					a = AutoModel(env, alnfile=osp.join(self.alignments_dir, additional_modeller_template + '.ali'), knowns=pdb_chain, sequence=uniprot, assess_methods=(assess.DOPE, assess.GA341), root_name=additional_modeller_template)                             
					a.starting_model = 1
					a.ending_model = 1
					a.max_molpdf = 1e8 # default is 1e7, increase for some not so great alignments/PDB files
					a.make()
				except:
					print('ERROR!', uniprot, pdb_chain)
			else:
				print('Model has already been found for:', uniprot, pdb_chain + '...')

	def get_additional_models(self):
		self.get_auth_label_dict()
		self.get_label_auth_dict()
		self.convert_blast_alignments()
		self.create_uniprot_ali_files()
		self.build_additional_models()

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	parser.add_argument('-c', '--scratch_dir') #/home/username/scratch
	parser.add_argument('-p', '--pdb_download_dir') #/home/username/scratch/pdb_cif
	args = parser.parse_args()

	print('Getting additional MODELLER models...')
	c = AdditionalModels(args.script_dir, args.scratch_dir, args.pdb_download_dir)
	c.get_additional_models()


if __name__=='__main__':
	main()