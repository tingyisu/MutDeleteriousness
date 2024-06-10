'''
Processes binary interactions in the IntAct dataset of human protein-protein interactions
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from Bio import SeqIO
import pickle
import urllib
import pandas as pd
from io import StringIO
import csv
from simple_tools import pickle_load, pickle_dump


class Data: 
    def __init__(self, orig_data_dir, processed_data_dir):
        # files to process
        self.ref_ppi_file = osp.join(orig_data_dir, 'intact.txt')
        # attributes to load 
        self.swissprot_ids = pickle_load(osp.join(processed_data_dir, 'swissprot_ids_list.pickle')) # list of reviewed Uniprotkb proteins
        self.swissprot_gene = pickle_load(osp.join(processed_data_dir, 'uniprot_gene_name_dict.pickle')) # key = swissprot ID, value = list of genes (put them all in file? pick one?)
        self.swissprot_seq = pickle_load(osp.join(processed_data_dir, 'swissprot_seq_dict.pickle')) # key = uniprot id, value = protein sequence; only for uniprot proteins in unique_uniprot
        # attributes to save
        self.pairs = [] # list of tuples interacting uniprot proteins
        self.unique_uniprot = [] # list of unique Uniprot proteins in self.pairs (for BLAST)
        self.gene_dict = {} # key = uniprot id, value = gene name
        self.intact_id_dict = {} # key = uniprot id, value = IntAct ID
        # files to write to/pickle attributes to
        interactome_data_dir = osp.join(processed_data_dir, 'interactome')
        self.intact_physical_interactions_file = osp.join(interactome_data_dir, 'intact_physical_interactions.tsv')
        self.intact_physical_interactions_processed_file = osp.join(interactome_data_dir, 'intact_physical_interactions_processed.tsv')
        self.pairs_file = osp.join(interactome_data_dir, 'intact_uniprot_pairs_physical.pickle')
        self.gene_dict_file = osp.join(interactome_data_dir, 'intact_gene_dict_physical.pickle')
        self.intact_id_dict_file = osp.join(interactome_data_dir, 'intact_id_dict_physical.pickle')
        self.uniprot_fasta_blast_file = osp.join(interactome_data_dir, 'intact_uniprot_sequences_blast.fasta')


    def get_physical_interactions(self):
        print('-----Getting physical interactions from intact.txt-----')
        # count number of interactions of a specific type
        physical_interactions = ['direct interaction', 'physical association', 'association']
        interaction_type_dict = {} # keys = interactions in physical_interactions, values = number of interactions
        num_physical_interactions = 0
        with open(self.ref_ppi_file, 'r') as f, open(self.intact_physical_interactions_file, 'w') as fout:
            header = f.readline()
            fout.write(header)
            for line in f:
                interaction_type = line.split('\t')[11]
                interaction = interaction_type.split('(')[1].split(')')[0]
                if interaction in physical_interactions:
                    num_physical_interactions += 1
                    fout.write(line)
                    if interaction not in interaction_type_dict:
                        interaction_type_dict[interaction] = 1
                    else:
                        interaction_type_dict[interaction] += 1
                    
        print('Total number of physical interactions:', num_physical_interactions)
        print('Number of direct interactions:', interaction_type_dict['direct interaction'])
        print('Number of physical associations:', interaction_type_dict['physical association'])
        print('Number of associations:', interaction_type_dict['association'])

    def process_intact_file(self):
        '''
        adapted from parse_IntAct_interactions in 
        https://github.com/MohamedGhadie/dispensable_ppi_content/blob/457bedf72716681a1074ac05cbca3eb18afb7599/code/text_tools.py#L17

        self.intact_physical_interactions_file contains interactions for all 16 species, so uses only swissprot reviewed human proteins to filter the interactions
        '''
        print('-----Finding pairs of interacting swissprot proteins in IntAct-----')

        with open(self.intact_physical_interactions_file, 'r') as f, open(self.intact_physical_interactions_processed_file, 'w') as fout:

            fout.write('\t'.join(["Gene1", 
                                "Gene2", 
                                "Protein1",
                                "Protein2",
                                "Protein1_IntAct_ID",
                                "Protein2_IntAct_ID",
                                "Interaction_ID"]) + '\n')
            i = c = 0
            next(f)
            for line in f:
                i += 1
                strsplit = line.split('\t')
                if strsplit[35].lower()=='false':
                    if (len(strsplit[0])>10) and (len(strsplit[1])>10):
                        if (strsplit[0][:10].lower()=='uniprotkb:') and (strsplit[1][:10].lower()=='uniprotkb:'):
                            protein1 = strsplit[0][10:]
                            protein2 = strsplit[1][10:]
                            if (protein1 in self.swissprot_ids) and (protein2 in self.swissprot_ids): # keep all interactions, including self interactions
                                IDsplit = list(map(str.strip, strsplit[2].split('|')))
                                hasIntactID = list(map(lambda x: x.find('intact:EBI-'), IDsplit))
                                try:
                                    ind = hasIntactID.index(0)
                                    intactID1 = IDsplit[ind][7:]
                                except ValueError:
                                    intactID1 = '-'
                                IDsplit = list(map(str.strip, strsplit[3].split('|')))
                                hasIntactID = list(map(lambda x: x.find('intact:EBI-'), IDsplit))
                                try:
                                    ind = hasIntactID.index(0)
                                    intactID2 = IDsplit[ind][7:]
                                except ValueError:
                                    intactID2 = '-'
                                IDsplit = list(map(str.strip, strsplit[13].split('|')))
                                hasIntactID = list(map(lambda x: x.find('intact:EBI-'), IDsplit))
                                try:
                                    ind = hasIntactID.index(0)
                                    interactionID = IDsplit[ind][7:]
                                except ValueError:
                                    interactionID = '-'
                                # get only human interactions using UniProtKB/SWISS-PROT (reviewed) proteins
                                gene1 = self.swissprot_gene[protein1] if protein1 in self.swissprot_gene else '-'
                                gene2 = self.swissprot_gene[protein2] if protein2 in self.swissprot_gene else '-'
                                # genes= [gene1, gene2]
                                c += 1
                                fout.write('\t'.join([','.join(gene1),
                                                    ','.join(gene2),
                                                    protein1,
                                                    protein2,
                                                    intactID1,
                                                    intactID2,
                                                    interactionID]) +  '\n')
                                # print('wrote:', genes, protein1, protein2)
        print('%d lines parsed from ' % i + 'IntAct reference file'
          ', %d PPIs extracted ' % c + 
          'and written to file ' + 'IntAct processed file')

    # gets self.pairs, self.gene_dict, self.intact_id_dict
    def get_pairs_gene_intact_id(self, incl_self, self_only):
        '''
        args:
            inc_self = T/F, T = include self-interactions
            self_only = T/F, T = only include self-interactions
        '''
        print('-----Finding pairs of Uniprot Interacting Proteins that appear more than once-----')
        pair_count_dict = {}
        with open(self.intact_physical_interactions_processed_file, 'r') as f:
            next(f) # skip header line
            for line in f:
                items = line.split('\t') # items[2] and items[3] are the proteins
                if items[2] not in self.gene_dict:
                    self.gene_dict[items[2]] = items[0]
                if items[3] not in self.gene_dict:
                    self.gene_dict[items[3]] = items[1]
                if items[2] not in self.intact_id_dict:
                    self.intact_id_dict[items[2]] = items[4]
                if items[3] not in self.intact_id_dict:
                    self.intact_id_dict[items[3]] = items[5]
                # count number of times pair occurs in self.intact_physical_interactions_processed_file
                if not incl_self and items[2] == items[3]: # to check whether to include self interactions
                    continue
                if self_only and items[2] != items[3]:
                    continue
                else:
                    if (items[2], items[3]) not in pair_count_dict and (items[3], items[2]) not in pair_count_dict:
                        pair_count_dict[(items[2], items[3])] = 1
                    elif (items[2], items[3]) in pair_count_dict:
                        pair_count_dict[(items[2], items[3])] += 1
                    else: # (items[3], items[2]) in pair_count_dict
                        pair_count_dict[(items[3], items[2])] += 1
        # check that there are no repeated pairs in reverse order i.e. (p1, p2) and (p2, p1) are both in pair_count_dict
        if not incl_self:
            for (p1, p2) in pair_count_dict:
                if (p2, p1) in pair_count_dict:
                    print('oh no! The reverse-order pair is included:', (p1, p2), (p2, p1))
        # only take interacting pairs that occur more than once in self.intact_physical_interactions_processed_file
        # 2 hits!
        for pair in pair_count_dict:
            if pair_count_dict[pair] > 1:
                self.pairs.append(pair)

    # --- the following functions get sequences for unique uniprot proteins for BLAST ---
    # gets unique uniprot proteins (for blast)
    def get_unique_uniprot(self):
        print('-----Getting list of unique uniprot IDs for BLAST-----')
        for (p1, p2) in self.pairs:
            self.unique_uniprot.extend([p1, p2])
        self.unique_uniprot = list(set(self.unique_uniprot))

    # create UniProt fasta file to run BLASTP
    def get_uniprot_fasta_blast(self):
        print('-----Writing Uniprot Seq to fasta BLAST format-----')
        with open(self.uniprot_fasta_blast_file, 'w') as f:
            for p in self.unique_uniprot:
                f.write('>' + str(p) + '\n')
                f.write(str(self.swissprot_seq[p]) + '\n')

    def process_intact_data(self):
        self.get_physical_interactions()
        self.process_intact_file() # get only human swissprot reviewed interacting protein pairs
        self.get_pairs_gene_intact_id(False, False) # remove self-interactions
        self.get_unique_uniprot()
        self.get_uniprot_fasta_blast()

        # pickle attributes to files
        pickle_dump(self.pairs, self.pairs_file)
        pickle_dump(self.gene_dict, self.gene_dict_file)
        pickle_dump(self.intact_id_dict, self.intact_id_dict_file)

        print('Number of unique Uniprot pairs found:', len(self.pairs))
        print('Number of unique Uniprot proteins in the interacting pairs:', len(self.unique_uniprot))


def main():
    script_dir = osp.dirname(__file__)
    orig_data_dir = osp.join(script_dir, '..', 'data', 'original')
    processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')

    d = Data(orig_data_dir, processed_data_dir)
    d.process_intact_data()
    

if __name__=='__main__':
    main()