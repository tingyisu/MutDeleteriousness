'''
Contains scripts to find interfacial residues and memoize them (reduces runtime for repeated PDB chain-pairs)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
import pickle
import os.path as osp
import argparse 
import os
import timeit

class Residues:
    def __init__(self, name):
        self.name = name
        self.memo_residues = {}
        self.memo_keys = {}
        self.pos_hom = {}
        self.pdb_chain_dict = {}
        self.num_chains = 0

    # read in memo_residues file
    def get_memo_from_file(self, memo_residues_file):
        if osp.exists(memo_residues_file):
            print('-----Loading residues from memo file-----')
            with open(memo_residues_file, 'r') as fname:
                for line in fname:
                    items = line.strip().split('\t')
                    struct_name, chain1_name, chain2_name = items[0], items[1], items[2]
                    if len(items) == 3: # no interfacial residues
                        self.memo_residues[(struct_name, chain1_name, chain2_name)] = [[],[]]
                    else:
                        residue1_str, residue2_str = items[3], items[4]
                        self.memo_residues[(struct_name, chain1_name, chain2_name)] = [residue1_str.split(','), residue2_str.split(',')]
        else:
            print('No memo residues file exists yet.')

    # get unique PDB chains
    # returns dict of key = PDB structure, value = dict of possibly interacting chains
    def get_pdb_chain_dict(self, pos_hom_file, inc_self, self_only, f_exists, pdb_chain_dict_file):
        '''
        args:
            inc_self = T/F, T = include self-interactions
            self_only = T/F, T = only include self-interactions
            f_exists = T/F, whether pickled pdb_chain_dict exists already
            pdb_chain_dict_file = name of file containing pickled self.pdb_chain_dict attribute
        '''
        print('-----Getting unique PDB chains-----')
        if f_exists:
            with open(pdb_chain_dict_file, 'rb') as f:
                self.pdb_chain_dict = pickle.load(f)
            self.num_chains = 0
            for key in self.pdb_chain_dict:
                self.num_chains += len(self.pdb_chain_dict[key])
            print('Total number of chains: ', self.num_chains)
        else:
            # read in self.pos_hom from file
            with open(pos_hom_file, 'rb') as fname_pos_hom:
                self.pos_hom = pickle.load(fname_pos_hom)

            num_uniprot_pairs = 0
            # iterate through self.pos_hom to get unique interacting chains
            for pair in self.pos_hom:

                p1, p2 = pair

                # if don't want to include self-interactions
                if not inc_self and p1 == p2:
                    continue

                # if want self interactions only
                if self_only and p1 != p2:
                    continue

                num_uniprot_pairs += 1

                for tup in self.pos_hom[pair]:
                    # get struct_name and chain names
                    struct_name = tup[0].split('_')[0]
                    chain1_name, chain2_name = tup[0].split('_')[1], tup[1].split('_')[1]
                    if chain1_name == chain2_name: print('same chain:', struct_name, chain1_name, chain2_name)
                    self.pdb_chain_dict.setdefault(struct_name, []).append((chain1_name, chain2_name))

            # get total number of unique chains as well as unique chains dict
            for struct_name in self.pdb_chain_dict:
                chain_list = list(set(self.pdb_chain_dict[struct_name]))
                # get unique chain_list
                new_chain_list = []
                for (chain1_name, chain2_name) in chain_list:
                    if (chain1_name, chain2_name) not in new_chain_list and (chain2_name, chain1_name) not in new_chain_list:
                        new_chain_list.append((chain1_name, chain2_name))
                self.pdb_chain_dict[struct_name] = new_chain_list
                self.num_chains += len(self.pdb_chain_dict[struct_name])

            print('Number of uniprot pairs that share at least one PDB structure:', num_uniprot_pairs)
            print('Number of unique interacting PDB chains:', self.num_chains)
            print('Number of PDB structures:', len(self.pdb_chain_dict))

    # do this only for IntAct (to avoid finding interfacial residues for the PDB chains-pairs that are also present in HI-union)
    def update_pdb_chain_dict(self, other_pdb_chain_dict_file):
        other_pdb_chain_dict = {}
        with open(other_pdb_chain_dict_file, 'rb') as f:
            other_pdb_chain_dict = pickle.load(f)

        # only keep interacting chains that aren't in other_pdb_chain_dict_file
        to_delete = [] # if have no chains leftover
        for struct in self.pdb_chain_dict:
            if struct in other_pdb_chain_dict:
                new_chain_list = []
                for (chain1_name, chain2_name) in self.pdb_chain_dict[struct]:
                    if (chain1_name, chain2_name) not in other_pdb_chain_dict[struct] and (chain2_name, chain1_name) not in other_pdb_chain_dict[struct]:
                        new_chain_list.append((chain1_name, chain2_name))
                if len(new_chain_list) == 0:
                    to_delete.append(struct)
                    # print('none left for:', struct)
                elif new_chain_list == self.pdb_chain_dict[struct]:
                    print('still no change:', struct)
                    # continue
                else:
                    self.pdb_chain_dict[struct] = new_chain_list

        # delete struct keys that have no chains left over
        for struct in to_delete:
            del self.pdb_chain_dict[struct]

        # update self.num_chains
        self.num_chains = 0
        for struct in self.pdb_chain_dict:
            self.num_chains += len(self.pdb_chain_dict[struct])

        print('Number of non-overlapping chains:', self.num_chains)

    def save_pdb_chain_dict(self, pdb_chain_dict_file):
        with open(pdb_chain_dict_file, 'wb') as f:
            pickle.dump(self.pdb_chain_dict, f) 

     # gets PDB structures and downloads them
    def download_pdb_structures(self, pdb_download_path):
        print('-----Downloading PDB structures-----')
        pdb = PDBList()
        # create pdb directory if it doesn't exist already
        if not osp.exists(pdb_download_path):
            os.mkdir(pdb_download_path)
        # download pdb structure if it doesn't exist in pdb directory already
        for struct in list(self.pdb_chain_dict.keys()):
            if not osp.exists(osp.join(pdb_download_path, struct + '.cif')):
                print(f'--Downloading {struct}--')
                pdb.retrieve_pdb_file(struct, pdir=pdb_download_path)

    def get_new_pdb_chain_dict(self):
        '''
        returns new dictionary that hasn't been looked at yet
        '''
        print('-----Getting new pdb chain dict-----')
        new_pdb_chain_dict = {}
        # get self.memo_keys
        for key in self.memo_residues:
            struct_name, chain1_name, chain2_name = key
            self.memo_keys.setdefault(struct_name, []).append((chain1_name, chain2_name))
        # update new_pdb_chain_dict with unseen residues
        self.num_chains = 0
        for key in self.pdb_chain_dict:
            if key not in self.memo_keys:
                new_pdb_chain_dict[key] = self.pdb_chain_dict[key]
                self.num_chains += len(self.pdb_chain_dict[key])
            else:
                unseen_chains = []
                for (chain1_name, chain2_name) in self.pdb_chain_dict[key]:
                    if (chain1_name, chain2_name) not in self.memo_keys[key] and (chain2_name, chain1_name) not in self.memo_keys[key]:
                        unseen_chains.append((chain1_name, chain2_name))
                if unseen_chains != []:
                    self.num_chains += len(unseen_chains)
                    new_pdb_chain_dict[key] = unseen_chains
        print('Number of unseen chains:', self.num_chains)
        # print(new_pdb_chain_dict['4awl'])
        return new_pdb_chain_dict, self.num_chains

    def get_residues(self, pdb_download_path, memo_residues_file, rerun):
        '''
        rerun = T/F, whether or not running again on compute canada
        '''
        print('------Finding residues-----')
        chain_dict = {}
        i = 0
        if not rerun:
            chain_dict = self.pdb_chain_dict
        else:
            chain_dict, _ = self.get_new_pdb_chain_dict()
        # return
        # create parser
        parser = MMCIFParser(QUIET=True) # suppress warnings
        pdb = PDBList()
        print('-----Updating memo_residues-----')
        with open(memo_residues_file, 'a') as fname_memo:
            for struct_name in chain_dict: # struct_name = pdb structure name
                print('Looking at structure:', struct_name)
                # get struct_name and chain names
                fname = struct_name + '.cif'
                # add just in case compute canada deletes .cif files
                structure = ""
                try:
                    structure = parser.get_structure(struct_name, osp.join(pdb_download_path, fname))
                except Exception as e: # redownload .cif file if need be
                    print(f'--Downloading {fname}--')
                    pdb.retrieve_pdb_file(struct_name, pdir=pdb_download_path)
                    structure = parser.get_structure(struct_name, osp.join(pdb_download_path, fname))
                # if not osp.exists(osp.join(pdb_download_path, fname)):
                #     print(f'--Downloading {fname}--')
                #     pdb.retrieve_pdb_file(struct_name, pdir=pdb_download_path)
                # get structure, model, and chains
                # structure = parser.get_structure(struct_name, osp.join(pdb_download_path, fname))
                model = structure[0]
                for chain_tup in chain_dict[struct_name]:
                    print('Chain:', struct_name, chain_tup)
                    chain1_name, chain2_name = chain_tup
                    i += 1
                    if (struct_name, chain1_name, chain2_name) not in self.memo_residues and (struct_name, chain2_name, chain1_name) not in self.memo_residues:
                        num_atoms = 0
                        start = timeit.default_timer() 
                        # get interfacial residues (residues that are within 5 angstroms of each other)
                        chain1, chain2 = model[chain1_name], model[chain2_name]
                        residue1_list, residue2_list = [], []
                        for residue1 in chain1:    
                            if residue1.get_full_id()[3][0] == " ": #ignore the heteroatoms 
                                res1_num = str(residue1.get_full_id()[3][1])
                                res1_id = str(residue1.get_full_id()[3][2]) # add in the auth_seq_id identifier (if exists), otherwise could have multiple label_seq_id for the same auth_seq_id (e.g. for auth_seq_id 156 on chain A of 1kmc --> has 156 and 156A)
                                res1 = res1_num
                                if res1_id != ' ':
                                    res1 = res1_num + res1_id
                                for residue2 in chain2:
                                    if residue2.get_full_id()[3][0] == " ": #ignore the heteroatoms 
                                        res2_num = str(residue2.get_full_id()[3][1])
                                        res2_id = str(residue2.get_full_id()[3][2]) # add in the auth_seq_id identifier, otherwise could have multiple label_seq_id for the same auth_seq_id (e.g. for auth_seq_id 156 on chain A of 1kmc --> has 156 and 156A)
                                        res2 = res2_num
                                        if res2_id != ' ':
                                            res2 = res2_num + res2_id
                                        # compute distance between all atoms
                                        try:
                                            break_loop = False
                                            # print(f'Looking at {res1_num}, {res2_num}')
                                            for atom1 in residue1: 
                                                for atom2 in residue2:
                                                    # ways to calculate euclidean distance
                                                    # 1. np.linalg.norm(residue1[atom1.name].get_coord() - residue2[atom2.name].get_coord())
                                                    # 2. residue1[atom1.name] - residue2[atom2.name] 
                                                    # atom1 - atom2 <= 5 --> might be faster
                                                    # check distances between x, y, z coordinates
                                                    # if atom1[0] - atom2[0] > 5 
                                                    num_atoms += 1
                                                    atom1_coords, atom2_coords = atom1.get_coord(), atom2.get_coord()
                                                    # don't keep searching if distance between each coordinate (x, y, z) > 5
                                                    if abs(atom1_coords[0] - atom2_coords[0]) > 5 or abs(atom1_coords[1] - atom2_coords[1]) > 5 or abs(atom1_coords[2] - atom2_coords[2]) > 5:
                                                        continue
                                                    else:
                                                        if atom1 - atom2 <= 5: # if the euclidean distance between any atoms in residues 1 & 2 <= 5, then they are interacting
                                                            if res1 not in residue1_list: residue1_list.append(res1)
                                                            if res2 not in residue2_list: residue2_list.append(res2)
                                                            break_loop = True
                                                            break
                                                if break_loop: break
                                        except KeyError:
                                            continue
                                                # print(residue1, residue2, distance)

                        # dynamic programming: save interfacial residues that have been looked at already
                        # write to memo file regardless of whether residue lists are both empty

                        stop = timeit.default_timer()
                        self.memo_residues[(struct_name, chain1_name, chain2_name)] = [residue1_list, residue2_list]
                        fname_memo.write('\t'.join([struct_name, chain1_name, chain2_name, ",".join(residue1_list), ",".join(residue2_list)]) + '\n')
                        fname_memo.flush() # need to flush so that will appear in output file immediately
                        print(str(i) + '/' + str(self.num_chains) + ' ' + '\t'.join([struct_name, chain1_name, chain2_name]) + '; time taken: ' + str(stop-start) + '; num atoms:' + str(num_atoms))

                    else:
                        print('ALREADY FOUND:' + str(i) + '/' + str(self.num_chains) + ' ' + '\t'.join([struct_name, chain1_name, chain2_name]))
