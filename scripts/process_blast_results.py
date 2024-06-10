'''
Processes BLASTP sequence alignment results
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
import pickle
import os.path as osp
import argparse 
import os

class Interactions:
    def __init__(self, name):
        self.name = name
        self.dict = {} # dict of blast results for interacting proteins, key = uniprot ID, value = list of PDB chains
        self.uniprot_pairs = [] # list of uniprot interacting pairs
        self.pos_hom = {} # possible homologs (in the same PDB structure), key = uniprot pair, value = list of tuples of PDB chains
    
    # creates dict of blast results for interacting proteins
    def create_dict(self, blast_file):
        print('-----Creating dict of blast results-----')
        f = open(blast_file, 'r')
        for line in f.readlines():
            items = line.strip().split('\t')
            self.dict.setdefault(items[0], []).append(items[1]) # add to dict

    # get uniprot interacting pairs from pickle file
    def get_unique_uniprot_pairs(self, pairs_pickle_file):
        print('-----Getting pairs of unique interacting Uniprot IDs-----')
        with open(pairs_pickle_file, 'rb') as fname:
            self.uniprot_pairs = pickle.load(fname)

    # finds matching PDB structures for a pair of 'interacting' proteins (uniprot ID)
    def match_PDB(self):
        print('-----Matching PDB structures for pairs of interacting proteins-----')
        for pair in self.uniprot_pairs:
            p1, p2 = pair
            if p1 in self.dict and p2 in self.dict:
                for pdb1 in self.dict[p1]: 
                    found = [(pdb1, pdb2) for pdb2 in self.dict[p2] if pdb1.split('_')[0] in pdb2 and pdb1 != pdb2]
                    if found != []: 
                        self.pos_hom.setdefault(pair, []).extend(found) # add to dict if the pair shares any of the same pdb structures      

    # finds the number of self-interactions in self.pos_hom
    def find_num_interactions(self):
        num_self = 0
        num_non_self = 0
        for (p1, p2) in self.pos_hom:
            if p1 == p2: 
                num_self += 1
            else: 
                num_non_self += 1
        print('Number of self interactions:', num_self)
        print('Number of non-self interactions:', num_non_self)
    
    # saves self.pos_hom dict to .pickle file
    def save_pos_hom(self, pos_hom_file):
        with open(pos_hom_file, 'wb') as fname:
            pickle.dump(self.pos_hom, fname)

                
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_directory') # e.g. /home/username/scratch/interfacial_residues/files
    args = parser.parse_args()

    print('-----HI-union-----')
    hiunion_blast_results_file = osp.join(args.output_directory, 'hiunion_blast_for_finding_interactions.tsv')
    hiunion_pairs_file = osp.join(args.output_directory, 'HI-union_unique_uniprot_pairs.pickle')
    hiunion_pos_hom_file = osp.join(args.output_directory, 'hiunion_pos_hom.pickle')

    hiunion = Interactions('hiunion')
    hiunion.create_dict(hiunion_blast_results_file)
    hiunion.get_unique_uniprot_pairs(hiunion_pairs_file)    
    hiunion.match_PDB()
    hiunion.save_pos_hom(hiunion_pos_hom_file)
    hiunion.find_num_interactions()

    print('-----IntAct-----')
    intact_blast_results_file = osp.join(args.output_directory, 'intact_blast_for_finding_interactions.tsv')
    intact_pairs_file = osp.join(args.output_directory, 'intact_uniprot_pairs_physical.pickle')
    intact_pos_hom_file = osp.join(args.output_directory, 'intact_pos_hom_physical.pickle')

    intact = Interactions('intact')
    intact.create_dict(intact_blast_results_file)
    intact.get_unique_uniprot_pairs(intact_pairs_file)
    intact.match_PDB()
    intact.save_pos_hom(intact_pos_hom_file)
    intact.find_num_interactions()


if __name__=='__main__':
    main()