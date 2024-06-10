'''
Script to run MODELLER on all alignments within a certain job
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os
import os.path as osp
import argparse
import pickle
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from simple_tools import pickle_load, check_create_dir, check_exist, pickle_dump

class Modeller:
	def __init__(self, script_dir, n):
		self.data_dir = osp.join(script_dir, 'files')
		self.alignments_dir = osp.join(script_dir, 'alignments')
		check_create_dir(self.alignments_dir)
		self.results_dir = osp.join(script_dir, 'results')
		check_create_dir(self.results_dir)
		self.blast_alignments = pickle_load(osp.join(self.alignments_dir, 'blast_alignments_for_modeller_' + n + '.pickle')) # list of (protein, pdb_chain) pairs
		self.discrepancies_between_seq_res_mmcif = []
		self.discrepancies_to_remove_pickle_file = osp.join(self.data_dir, 'discrepancies_to_remove_' + n + '.pickle')

	def build_models(self):
		for (uniprot, pdb_chain) in self.blast_alignments:
			print('Finding model for:', uniprot, pdb_chain + '...')
			# run if haven't found model already
			fname = '_'.join([uniprot, pdb_chain])
			model_name = osp.join(self.results_dir, fname + '.B99990001.pdb')
			if not check_exist(model_name):
				try:
					fname = '_'.join([uniprot, pdb_chain])
					log.minimal()
					env = Environ()
					env.io.hetatm = False # (NOTE: MODELLER seems to take HETATM MSE regardless of whether this is set to True/False)! It converts MSE to MET! https://salilab.org/archives/modeller_usage/2009/msg00185.html
					a = AutoModel(env, alnfile=osp.join(self.alignments_dir, fname + '.ali'), knowns=pdb_chain, sequence=uniprot, assess_methods=(assess.DOPE, assess.GA341), root_name=fname)                             
					a.starting_model = 1
					a.ending_model = 1
					a.max_molpdf = 1e8 # default is 1e7, increase for some not so great alignments/PDB files
					a.make()
				except:
					print('ERROR! Discrepancies between PDB SeqRes and mmCIF for', uniprot, pdb_chain)
					self.discrepancies_between_seq_res_mmcif.append((uniprot, pdb_chain))
			else:
				print('Model has already been found for:', uniprot, pdb_chain + '...')
		pickle_dump(self.discrepancies_between_seq_res_mmcif, self.discrepancies_to_remove_pickle_file)

	# def model_evaluation(self):



def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-yxia/username/modeller
	parser.add_argument('-n', '--num_part') # integer
	args = parser.parse_args()

	m = Modeller(args.script_dir, args.num_part)
	m.build_models()

if __name__=='__main__':
	main()