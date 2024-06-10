'''
Simple functions that are used repeatedly
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import os.path as osp
import os
from Bio.PDB import PDBList

def pickle_load(pickle_file):
	with open(pickle_file, 'rb') as f:
		contents = pickle.load(f)
		return contents

def pickle_dump(contents, pickle_file):
	with open(pickle_file, 'wb') as f:
		pickle.dump(contents, f)

# check whether a file/directory exists or not
def check_exist(dirname):
	return osp.exists(dirname)

# creates a directory if it doesn't exist already
def check_create_dir(dirname):
	if not osp.exists(dirname):
		os.mkdir(dirname)

# remove directory named dirname if it exists
def check_remove(dirname):
	if osp.exists(dirname):
		os.remove(dirname)

# writes a list of lists to a .tsv file 
# items in the list are separated by \n
# items in the sublists are separated by \t
def write_list_of_lists_to_tsv(contents, file):
	with open(file, 'w') as f:
		for line in contents:
			line = [str(item) for item in line]
			f.write('\t'.join(line) + '\n')

# reads in items from newline delimited file
# and returns a list of the items
def get_list_from_file(newline_delimited_file):
	items = []
	with open(newline_delimited_file, 'r') as f:
		items = [line.strip() for line in f]
	return items

# read in list with header from newline delimited file
def read_list_with_header_from_file(file):
	items = []
	with open(file, 'r') as f:
		next(f)
		items = [line.strip() for line in f]
	return items
		
def write_list_to_file(list_to_write, file):
	with open(file, 'w') as f:
		for item in list_to_write:
			f.write(item + '\n')

# update self.missense_info and self.header_dict
def get_missense_info(missense_file):
	'''
	returns:
	1) header_dict = {'column1_name': column1_index, 'column2_name': column2_index, ...}
	2) missense_info = [[line1_column1_entry, line1_column2_entry,...], [line2_column1_entry, line2_column2_entry,...]]
	'''
	header_dict = {}
	missense_info = []
	with open(missense_file, 'r') as f:
		column_names = f.readline().strip().split('\t')
		missense_info.append(column_names) # missense_info also contains a list of the header column names
		i = 0
		for name in column_names:
			header_dict[name] = i
			i += 1
		for line in f:
			missense_info.append(line.strip().split('\t'))
	return header_dict, missense_info

# write to file
def write_missense_info_to_file(missense_info, missense_file):
	with open(missense_file, 'w') as f:
		for line in missense_info:
			line = [str(item) for item in line]
			f.write('\t'.join(line) + '\n')

# download PDB structure in mmcif format
# gets PDB structure and downloads it
def download_pdb_structure(pdb_structure, pdb_download_path):
	print('-----Downloading PDB structures in mmCIF format-----')
	pdb = PDBList()
	# create pdb directory if it doesn't exist already
	check_create_dir(pdb_download_path)
	# download pdb structure if it doesn't exist in pdb directory already
	if not check_exist(osp.join(pdb_download_path, pdb_structure + '.cif')):
		pdb.retrieve_pdb_file(pdb_structure, pdir=pdb_download_path)

# download PDB structures in mmcif format
# gets PDB structures and downloads them
def download_pdb_structures(pdb_structures_list, pdb_download_path):
	print('-----Downloading PDB structures in mmCIF format-----')
	pdb = PDBList()
	# create pdb directory if it doesn't exist already
	check_create_dir(pdb_download_path)
	# download pdb structure if it doesn't exist in pdb directory already
	for struct in pdb_structures_list:
		if not check_exist(osp.join(pdb_download_path, struct + '.cif')):
			pdb.retrieve_pdb_file(struct, pdir=pdb_download_path)

# create dictionary with keys & lists as values
def initialize_dict_of_lists(keys):
	initialized_dict = {}
	for key in keys:
		initialized_dict[key] = []
	return initialized_dict

# combines multiple lists into one list (removes overlapping items in lists)
def combine_lists(list_of_lists):
	all_items = []
	for l in list_of_lists:
		all_items += l
	return list(set(all_items))