'''
List .ali files with residue position discrepancies (MODELLER failed to build models for these .ali files due to such discrepancies)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse
import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump

class Discrepancies:
	def __init__(self, script_dir):
		self.data_dir = osp.join(script_dir, 'files')
		self.discrepancies_to_fix = []
		self.discrepancies_to_fix_pickle_file = osp.join(self.data_dir, 'rerun_ali_that_had_discrepancies.pickle')
		self.discrepancies_ali_files = []
	
	def list_discrepancies(self):
		discrepancy_pickle_files = [osp.join(self.data_dir, f) for f in os.listdir(self.data_dir) if 'discrepancies_to_remove' in f]
		if discrepancy_pickle_files != []:
			# print(discrepancy_pickle_files)
			for file in discrepancy_pickle_files:
				discrepancy_list = pickle_load(file)
				self.discrepancies_to_fix.extend(discrepancy_list)
			if len(set(self.discrepancies_to_fix)) != len(self.discrepancies_to_fix):
				print('There were repeated models being run...')
			else:
				pickle_dump(self.discrepancies_to_fix, self.discrepancies_to_fix_pickle_file) # pickle to file (for run_modeller_uncompleted.py)
		
				# print the .ali files
				for pair in self.discrepancies_to_fix:
					uniprot, template = pair
					self.discrepancies_ali_files.append(uniprot + '_' + template + '.ali')
				print('Number of protein-template chain pairs that have discrepancies between their PDB SeqRes and mmCIF:', len(self.discrepancies_to_fix))
				print('Fix the following discrepancies in the .ali files:')
				print(self.discrepancies_ali_files)

		else:
			print('No protein-template chain pairs that have discrepancies between their PDB SeqRes and mmCIF.')
			print('No need to change anything manually. Safe to continue.')
	
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-s', '--script_dir')  #/home/username/projects/def-*/username/modeller
	args = parser.parse_args()

	d = Discrepancies(args.script_dir)
	d.list_discrepancies()

if __name__=='__main__':
	main()
