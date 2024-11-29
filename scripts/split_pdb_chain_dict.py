'''
Splits PDB chain-pairs into groups and pickles them to different files (so that these groups can be run in parallel)
The HI-union and IntAct data are split separately
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import argparse
import numpy as np
from memo_residues import Residues
import os.path as osp
from simple_tools import pickle_dump, pickle_load

def split_pdb_chain_dict(pdb_chain_dict_file, num, num_chains):
	pdb_chain_dict = pickle_load(pdb_chain_dict)

	# split
	fname = pdb_chain_dict_file.split('.')[0]
	num_per_part = int(num_chains / num) 
	i = 0
	curr_part_num = 0
	curr_part_dict = {}
	number_to_save = [num_per_part for i in range(num)]
	number_to_save[len(number_to_save)-1] = num_chains - (num_per_part * (num-1))
	print(number_to_save)
	for key in pdb_chain_dict:
		for chain_pair in pdb_chain_dict[key]:
			if i < number_to_save[curr_part_num]:
				curr_part_dict.setdefault(key, []).append(chain_pair)
			else:
				print('Finished saving part', str(curr_part_num+1), 'of', str(i), 'chain pairs')
				pickle_dump(curr_part_dict, fname + '_split' + str(curr_part_num+1) + '.pickle')
				# print(curr_part_dict.keys())
				i = 0
				curr_part_num += 1
				curr_part_dict = {}
				curr_part_dict.setdefault(key, []).append(chain_pair)
			i += 1
	# for the last iteration
	print('Finished saving part', str(curr_part_num+1), 'of', str(i), 'chain pairs')
	pickle_dump(curr_part_dict, fname + '_split' + str(curr_part_num+1) + '.pickle')


# for combining leftover pdb_chain_dicts that didn't finish running in their previous job
def combine_pdb_chain_dict(combined_pdb_chain_dict_file, pdb_chain_dict_file_list):
	num_chains = 0
	combined_pdb_chain_dict = {}
	for file in pdb_chain_dict_file_list:
		curr_dict = pickle_load(file)
		for key in curr_dict:
			num_chains += len(curr_dict[key])
			combined_pdb_chain_dict.setdefault(key, []).extend(curr_dict[key]) # add to dict
		print('addition of file:', file, num_chains)
	pickle_dump(combined_pdb_chain_dict, combined_pdb_chain_dict_file)

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-o', '--output_directory') # /home/username/scratch/interfacial_residues/files
	parser.add_argument('-u', '--num_splits_for_hiunion') # number of slurm jobs to split HI-union data into
	parser.add_argument('-i', '--num_splits_for_intact') # number of slurm jobs to split IntAct data into
	args = parser.parse_args()

	name = 'hiunion'
	hiunion_memo_residues_file = osp.join(args.output_directory, name + "_memo_residues.tsv")
	hiunion_pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict.pickle')
	hiunion_new_pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict_new.pickle') # updated after getting memo residues

	hiunion = Residues(name)
	hiunion.get_memo_from_file(hiunion_memo_residues_file)
	hiunion.get_pdb_chain_dict('', False, False, True, hiunion_pdb_chain_dict_file)
	hiunion_new_pdb_chain_dict, hiunion_num_chains = hiunion.get_new_pdb_chain_dict()
	pickle_dump(hiunion_new_pdb_chain_dict, hiunion_new_pdb_chain_dict_file)
	split_pdb_chain_dict(hiunion_new_pdb_chain_dict_file, int(args.num_splits_for_hiunion), hiunion_num_chains)

	name = 'intact'
	intact_memo_residues_file = osp.join(args.output_directory, name + "_memo_residues.tsv")
	intact_pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict.pickle')
	intact_new_pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict_new.pickle') # updated after getting memo residues

	intact = Residues(name)
	intact.get_memo_from_file(intact_memo_residues_file)
	intact.get_pdb_chain_dict('', False, False, True, intact_pdb_chain_dict_file)
	intact_new_pdb_chain_dict, intact_num_chains = intact.get_new_pdb_chain_dict()
	pickle_dump(intact_new_pdb_chain_dict, intact_new_pdb_chain_dict_file)
	split_pdb_chain_dict(intact_new_pdb_chain_dict_file, int(args.num_splits_for_intact), intact_num_chains)

if __name__=='__main__':
	main()

