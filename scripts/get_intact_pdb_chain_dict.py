'''
Downloads PDB structures for interacting pairs of PDB chains that are homologous to interacting protein pairs in IntAct
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

    name = 'intact'
    print('-----IntAct-----')
    memo_residues_file = osp.join(args.output_directory, name + "_memo_residues.tsv")
    pos_hom_file = osp.join(args.output_directory, name + '_pos_hom_physical.pickle')
    pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict.pickle')
    hiunion_pdb_chain_dict_file = osp.join(args.output_directory, 'hiunion_pdb_chain_dict.pickle')
    output_num_file = osp.join(args.output_directory, name + '_num_chain.txt')
    
    intact = Residues(name)
    intact.get_pdb_chain_dict(pos_hom_file, False, False, False, pdb_chain_dict_file)
    intact.update_pdb_chain_dict(hiunion_pdb_chain_dict_file) # only for IntAct (b/c can remove all overlapping pdb_chain tuples from HI-Union)
    intact.save_pdb_chain_dict(pdb_chain_dict_file)
    intact.download_pdb_structures(args.pdb_download_path)


if __name__=='__main__':
    main()