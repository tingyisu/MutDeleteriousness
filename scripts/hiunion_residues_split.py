'''
Finds interfacial residues between each pair of PDB chains within each HI-union SLURM job
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
    parser.add_argument('-o', '--output_directory') # e.g. /home/username/scratch/interfacial_residues
    parser.add_argument('-p', '--part_number') # job number (ranging from 0 to n, n being the number of splits the HI-union data was split into)
    args = parser.parse_args()

    name = 'hiunion'
    print('-----Finding interfacial residues for HI-Union part' + args.part_number + '-----')
    memo_residues_file = osp.join(args.output_directory, name + "_memo_residues_split" + args.part_number + '.tsv')
    pdb_chain_dict_file = osp.join(args.output_directory, name + '_pdb_chain_dict_new_split' + args.part_number + '.pickle')
    hiunion = Residues(name)
    hiunion.get_memo_from_file(memo_residues_file)
    hiunion.get_pdb_chain_dict('', False, False, True, pdb_chain_dict_file)
    hiunion.get_residues(args.pdb_download_path, memo_residues_file, True) # keep rerun to True on default (even if is False)


if __name__=='__main__':
    main()