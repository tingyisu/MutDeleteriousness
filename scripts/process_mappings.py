'''
Processes mappings necessary for processing & mapping mutations to the structural interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
from Bio import SeqIO
from simple_tools import pickle_dump, get_missense_info, pickle_load, get_list_from_file

'''
pickles the following mappings:
(2) RefSeq protein -> uniprot
(3) RefSeq mRNA -> RefSeq protein
(6) RefSeq protein -> seq
(8) (RefSeq prot, gene_name) -> UniProt

NOTES: 
finds UniProt proteins based on gene name and RefSeq protein ID
only takes canonical UniProt IDs (those without '-*' in their IDs) that are SwissProt-reviewed
'''

class Mappings:
	def __init__(self, data_dir, orig_data_dir):
		self.uniprot_mappings_file = osp.join(orig_data_dir, 'HUMAN_9606_idmapping.dat')
		self.lrg_ref_seq_gene_file = osp.join(orig_data_dir, 'LRG_RefSeqGene')
		self.ref_seq_protein_seq_file = osp.join(data_dir, 'merged_ref_seq_protein.fa')
		self.swissprot_fasta_file = osp.join(orig_data_dir, 'uniprot_reviewed_human_proteome.fasta')
		self.data_dir = data_dir
		# attributes to get
		self.swissprot_ids = get_list_from_file(osp.join(orig_data_dir, 'uniprot_reviewed_human_proteome.list')) # list of reviewed Uniprotkb proteins
		self.swissprot_seq_dict = {} # key = swissprot id, value = sequence
		self.refseq_prot_uniprot_dict = {} # key = RefSeq protein ID, value = list of UniProt protein ID(s)
		self.refseq_prot_seq_dict = {} # key = RefSeq protein ID, value = protein sequence
		self.refseq_mrna_refseq_prot_dict = {} # key = RefSeq mRNA ID, value = RefSeq protein ID
		self.uniprot_gene_name_dict = {} # key = UniProt protein ID, value = list of Gene names
		self.refseq_prot_gene_name_uniprot_dict = {} # key = (RefSeq protein ID, gene name), value = UniProt protein ID
		# files to pickle to
		self.swissprot_ids_list_pickle_file = osp.join(self.data_dir, 'swissprot_ids_list.pickle')
		self.swissprot_seq_dict_pickle_file = osp.join(self.data_dir, 'swissprot_seq_dict.pickle')
		self.refseq_prot_uniprot_dict_pickle_file = osp.join(self.data_dir, 'refseq_prot_uniprot_dict.pickle')
		self.refseq_prot_seq_dict_pickle_file = osp.join(self.data_dir, 'refseq_prot_seq_dict.pickle')
		self.refseq_mrna_refseq_prot_dict_pickle_file = osp.join(self.data_dir, 'refseq_mrna_refseq_prot_dict.pickle')
		# ---
		self.uniprot_gene_name_dict_pickle_file = osp.join(self.data_dir, 'uniprot_gene_name_dict.pickle')
		self.refseq_prot_gene_name_uniprot_dict_pickle_file = osp.join(self.data_dir, 'refseq_prot_gene_name_uniprot_dict.pickle')

	# updates self.swissprot_seq_dict
	# checks that self.swissprot_ids are all keys in self.swissprot_seq_dict
	# then pickles self.swissprot_seq_dict and  self.swissprot_ids to files (so that don't need to process again for IntAct) (no need for querying now)
	def get_swissprot_id_seq(self):
		print('------Getting swissprot sequences and pickling the ids list and seqs dict-----')
		with open(self.swissprot_fasta_file, 'r') as f:
			fasta_seqs = SeqIO.parse(f,'fasta')
			for fasta_seq in fasta_seqs:
				# print(fasta_seq.id.split('|')[1])
				swissprot = fasta_seq.id.split('|')[1]
				self.swissprot_seq_dict[swissprot] = str(fasta_seq.seq)

		# check keys in self.swssiprot_seq_dict
		# swissprot_ids = list(self.swissprot_seq_dict.keys())
		# print(sorted(swissprot_ids) == sorted(self.swissprot_ids))

		# dump to files
		pickle_dump(self.swissprot_ids, self.swissprot_ids_list_pickle_file)
		pickle_dump(self.swissprot_seq_dict, self.swissprot_seq_dict_pickle_file)

	# use to update self.refseq_prot_seq_dict
	# get seq_dicts first, so that can retrieve protein ID -> uniprot mappings for only those protein IDs w/ known protein sequences
	def get_protein_seq(self, name, file):
		prot_id_seq_dict = {}
		print('-----Getting protein sequences for', name,  'protein IDs-----')
		with open(file, 'r') as f:
			fasta_seqs = SeqIO.parse(f,'fasta')
			for fasta_seq in fasta_seqs:
				if fasta_seq.id in prot_id_seq_dict:
					if fasta_seq.seq != prot_id_seq_dict[fasta_seq.id]:
						print('Different protein sequence!!!!', fasta_seq.id)
				else:
					prot_id_seq_dict[fasta_seq.id] = fasta_seq.seq
		return prot_id_seq_dict


	# update self.refseq_mrna_refseq_prot_dict
	def get_refseq_mrna_refseq_prot(self):
		print('-----Getting RefSeq mRNA to protein ID mappings-----')
		header_dict, refseq_info = get_missense_info(self.lrg_ref_seq_gene_file)
		# print(len(refseq_info))

		for i in range(1, len(refseq_info)):
			mrna_acc = refseq_info[i][header_dict['RNA']]
			prot_acc = refseq_info[i][header_dict['Protein']] 
			# gene_name = items[header_dict['Symbol']]
			if prot_acc in self.refseq_prot_seq_dict: # take only the entries that have protein sequences
				# check that duplicates contain the same info (same prot_acc and gene_name)
				if mrna_acc in self.refseq_mrna_refseq_prot_dict:
					# print('Have already seen this mRNA accession!', line)
					if prot_acc != self.refseq_mrna_refseq_prot_dict[mrna_acc]:
						print('Different protein accession:', mrna_acc, prot_acc, self.refseq_mrna_refseq_prot_dict[mrna_acc])
				else:
					self.refseq_mrna_refseq_prot_dict[mrna_acc] = prot_acc

	# update self.refseq_prot_uniprot_dict 
	def get_uniprot_mappings(self):
		print('-----Getting UniProt mappings for RefSeq protein IDs-----')
		with open(self.uniprot_mappings_file, 'r') as f:
			for line in f:
				items = line.strip().split('\t')
				uniprot_protein = items[0]
				category = items[1]
				other_id = items[2]

				if '-' not in uniprot_protein and uniprot_protein in self.swissprot_ids: # take only canonical uniprotkb ID that is swissprot reviewed
					if category == 'RefSeq':
						if other_id in self.refseq_prot_seq_dict: # make sure that refseq_prot has an existing protein sequence (i.e. exists in self.refseq_prot_seq_dict)
							if other_id in self.refseq_prot_uniprot_dict: # already exists, duplicate entry
								self.refseq_prot_uniprot_dict[other_id].append(uniprot_protein)
								# if uniprot_protein != self.refseq_prot_uniprot_dict: # make sure uniprot proteins are the same in duplicate entries
								# 	print('Different uniprot ID for RefSeq prot:', other_id, self.refseq_prot_uniprot_dict[other_id], uniprot_protein)
							else:
								self.refseq_prot_uniprot_dict[other_id] = [uniprot_protein]
					elif category == 'Gene_Name':
						if uniprot_protein in self.uniprot_gene_name_dict:
							self.uniprot_gene_name_dict[uniprot_protein].append(other_id)
							# if prev_gene_name != other_id:
							# 	print('Different gene name for uniprot protein:', uniprot_protein, prev_gene_name, other_id)
						else:
							self.uniprot_gene_name_dict[uniprot_protein] = [other_id]
					else:
						continue
				else:
					continue

	# for updating self.refseq_prot_gene_name_uniprot_dict
	def create_tuple_dict(self, name, prot_uniprot_dict, uniprot_gene_name_dict):
		print('-----Getting tuple dict for', name + '-----')
		tuple_dict = {}
		for prot_id in prot_uniprot_dict:
			for uniprot_protein in prot_uniprot_dict[prot_id]:
				if uniprot_protein in uniprot_gene_name_dict:
					for gene_name in uniprot_gene_name_dict[uniprot_protein]:
						tuple_dict[(prot_id, gene_name)] = uniprot_protein
				else:
					continue
		return tuple_dict



	def get_and_pickle_all_mappings(self):
		# get seq_dicts first, so that can retrieve protein ID -> uniprot mappings for only those protein IDs w/ known protein sequences
		self.refseq_prot_seq_dict = self.get_protein_seq('RefSeq', self.ref_seq_protein_seq_file)

		print('Number of RefSeq protein sequences:', len(self.refseq_prot_seq_dict))

		pickle_dump(self.refseq_prot_seq_dict, self.refseq_prot_seq_dict_pickle_file)

		self.get_uniprot_mappings()

		# check
		print('Number of RefSeq/UniProt proteins in self.refseq_prot_uniprot_dict:', len(self.refseq_prot_uniprot_dict))
		print('Number of UniProt protein IDs with gene names:', len(self.uniprot_gene_name_dict))

		pickle_dump(self.refseq_prot_uniprot_dict, self.refseq_prot_uniprot_dict_pickle_file)
		pickle_dump(self.uniprot_gene_name_dict, self.uniprot_gene_name_dict_pickle_file)


		self.refseq_prot_gene_name_uniprot_dict = self.create_tuple_dict('RefSeq', self.refseq_prot_uniprot_dict, self.uniprot_gene_name_dict)


		# check
		print('Number of RefSeq + gene_name tuples in self.refseq_prot_gene_name_uniprot_dict:', len(self.refseq_prot_gene_name_uniprot_dict))

		pickle_dump(self.refseq_prot_gene_name_uniprot_dict, self.refseq_prot_gene_name_uniprot_dict_pickle_file)

		self.get_refseq_mrna_refseq_prot()

		print('Number of RefSeq mRNA to protein mappings:', len(self.refseq_mrna_refseq_prot_dict))
		pickle_dump(self.refseq_mrna_refseq_prot_dict, self.refseq_mrna_refseq_prot_dict_pickle_file)

		self.get_swissprot_id_seq()


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')
	orig_data_dir = osp.join(script_dir, '..', 'data', 'original')

	m = Mappings(data_dir, orig_data_dir)
	m.get_and_pickle_all_mappings()

if __name__=='__main__':
	main()
