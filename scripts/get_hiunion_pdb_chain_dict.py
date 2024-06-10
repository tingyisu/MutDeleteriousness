'''
Downloads PDB structures for interacting pairs of PDB chains that are homologous to interacting protein pairs in HI-union
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
import pickle
import os.path as osp
import argparse 
import os
from memo_residues import Residues

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--pdb_download_path') # e.g. /home/username/scratch/pdb_cif
    parser.add_argument('-o', '--output_directory') # e.g. /home/username/scratch/interfacial_residues/files
    args = parser.parse_args()

    name = 'hiunion'
    print('-----HI-union-----')
    memo_residues_file = osp.join(args.output_directory, name + "_memo_residues.tsv")
    pos_hom_file = osp.join(args.output_directory, name + '_pos_hom.pickle')
    pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict.pickle')
    output_num_file = osp.join(args.output_directory, name + '_num_chain.txt')

    hiunion = Residues(name)
    hiunion.get_pdb_chain_dict(pos_hom_file, False, False, False, pdb_chain_dict_file)
    hiunion.save_pdb_chain_dict(pdb_chain_dict_file)
    hiunion.download_pdb_structures(args.pdb_download_path)    

if __name__=='__main__':
    main()