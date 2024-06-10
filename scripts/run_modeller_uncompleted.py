'''
Rerun MODELLER on alignments that failed to produce models
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse
import os.path as osp
import os
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
from simple_tools import pickle_load, check_exist

class Uncompleted:
	def __init__(self, script_dir):
		self.data_dir = osp.join(script_dir, 'files')
		self.alignments_dir = osp.join(script_dir, 'alignments')
		self.results_dir = osp.join(script_dir, 'results')
		discrepancies_pickle_file = osp.join(self.data_dir, 'rerun_ali_that_had_discrepancies.pickle')
		self.discrepancies_to_remove = []
		if check_exist(discrepancies_pickle_file):
			self.discrepancies_to_remove = pickle_load(discrepancies_pickle_file) # if there are alignments w/ discrepancies between PDB mmCIF & PDB SeqRes
	
	# get alignments that failed to produce structures b/c they exceeded max_molpdf or the molpdf was nan
	# rerun using run_modeller_uncompleted.py with a change in the modeller optimization method
	# also rerun alignments that failed due to discrepancies between PDB SeqRes & mmCIF
	def run_uncompleted_modeller_structures(self): 
		ali_pairs = [f.split('.')[0] for f in os.listdir(self.alignments_dir) if '.ali' in f]
		modeller_pdb_pairs = [f.split('.')[0] for f in os.listdir(self.results_dir) if '.pdb' in f]

		print('Number of alignment pairs:', len(ali_pairs))
		print('Number of modeller PDB pairs:', len(modeller_pdb_pairs))

		uncompleted = set(ali_pairs) - set(modeller_pdb_pairs)
		print('Length uncompleted pairs:', len(uncompleted))

		discrepancies_list = ['_'.join([uniprot, pdb_chain]) for (uniprot, pdb_chain) in self.discrepancies_to_remove]

		# run for models with initial alignment discrepancies (that have since been manually fixed)
		for (uniprot, pdb_chain) in self.discrepancies_to_remove:
			self.build_model(uniprot, pdb_chain)

		# run for models that have molpdf = nan or > 1e8
		to_run = list(uncompleted - set(discrepancies_list))

		if to_run != []:
			print(to_run)
			for pair in to_run:
				uniprot, pdb_struct, chain = pair.split('_')
				pdb_chain = '_'.join([pdb_struct, chain])
				self.build_model_change_optimization(uniprot, pdb_chain)

	# run this for models with initial alignment discrepancies
	def build_model(self, uniprot, pdb_chain):
		print('Finding model for:', uniprot, pdb_chain + '...')
		# run if haven't found model already
		fname = '_'.join([uniprot, pdb_chain])
		model_name = osp.join(self.results_dir, fname + '.B99990001.pdb')
		# if not check_exist(model_name):
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
			print('ERROR!', uniprot, pdb_chain)


	# try this for models that have molpdf = nan or > 1e8
	# or models with ****** in some of the atom coordinates
	def build_model_change_optimization(self, uniprot, pdb_chain): # use for uncompleted run (try this)
		# run if haven't found model already
		fname = '_'.join([uniprot, pdb_chain])
		model_name = osp.join(self.results_dir, fname + '.B99990001.pdb')
		# if not check_exist(model_name):
		try:
			fname = '_'.join([uniprot, pdb_chain])
			log.minimal()
			env = Environ()
			env.io.hetatm = False # (MODELLER seems to take HETATM MSE regardless of whether this is set to True/False)! It converts MSE to MET!
			a = AutoModel(env, alnfile=osp.join(self.alignments_dir, fname + '.ali'), knowns=pdb_chain, sequence=uniprot, assess_methods=(assess.DOPE, assess.GA341), root_name=fname)                             
			a.starting_model = 1
			a.ending_model = 1

			# change default optimization method is cannot get good model under max_molpdf

			# Very thorough VTFM optimization:
			a.library_schedule = autosched.slow
			a.max_var_iterations = 300
			# Thorough MD optimization:
			a.md_level = refine.slow
			# Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
			a.repeat_optimization = 2

			print('max_molpdf is:', a.max_molpdf) # can also try increasing a.max_molpdf if it still doesn't work
			a.make()
		except:
			print('ERROR!')


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	args = parser.parse_args()

	u = Uncompleted(args.script_dir)
	u.run_uncompleted_modeller_structures()


if __name__=='__main__':
	main()
