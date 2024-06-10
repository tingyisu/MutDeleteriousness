'''
Creates final interactions file (keep interactions with >=50% of interfacial residues mapping to their corresponding SwissProt/Uniprot proteins)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import argparse
import os.path as osp
from Bio.PDB import PDBList
import os
from simple_tools import pickle_load, pickle_dump

'''
* checks whether >= 50% of interfacial residues of one pdb chain map to its aligned uniprot protein
* first takes interfacial residues of PDB templates in auth format and converts them to label format
* then iterates through the BLAST alignment of the pdb template and corresponding uniprot protein, and checks whether at the 
* interfacial residue position has aligned amino acids on both the PDB template and uniprot protein
* then saves this interfacial residues based on the uniprot alignment position and checks whether >= 50% of all interfacial residues are mapped
* the mapped interfacial residues in UniProt alignment position are saved to *_final_interfacial_residues.tsv
* need to save the following mappings PDB label -> UniProt position
'''

class StructuralInteractome:
	def __init__(self, name, data_dir):
		self.name = name
		self.auth_label_dict = pickle_load(osp.join(data_dir, self.name + "_auth_label_dict.pickle")) # dict of dict of dict of key = auth_seq_id, value = label_seq_id, (PDB - pdb_seq_res.txt and biopython both use auth as convention)
		self.blast_file = osp.join(data_dir, self.name + "_blast_results_with_seq_eval-5.tsv")
		if self.name == 'intact':
			self.interactions_file = osp.join(data_dir, self.name + "_interfacial_residues_physical.tsv")
			self.final_interactions_file = osp.join(data_dir, self.name + "_final_interfacial_residues_physical.tsv")
			self.best_blast_alignment_file = osp.join(data_dir, self.name + "_blast_best_alignments_physical.tsv")
			self.best_blast_alignment_pickle_file = osp.join(data_dir, self.name + "_blast_best_alignments_physical.pickle")
		else:
			self.interactions_file = osp.join(data_dir, self.name + "_interfacial_residues.tsv")
			self.final_interactions_file = osp.join(data_dir, self.name + "_final_interfacial_residues.tsv")
			self.best_blast_alignment_file = osp.join(data_dir, self.name + "_blast_best_alignments.tsv")
			self.best_blast_alignment_pickle_file = osp.join(data_dir, self.name + "_blast_best_alignments.pickle")
		self.uniprot_pdb_chains = pickle_load(osp.join(data_dir, self.name + '_uniprot_pdb_chains.pickle'))
		self.blast_dict = {} # key = (uniprot, pdb chain), value = list of list of blast results
		self.best_blast_alignment = {} # key = (uniprot, pdb chain), value = list of list of best blast alignment result(s) (picked based on % mapped)
		# to pickle
		self.label_uniprot_dict = {} # dict of label-> uniprot mappings of interfacial residues that are mapped between PDB template and UniProt protein, key = ((uniprot1, pdb_template1), (uniprot2, pdb_template2)), value = [{label_seq_id1: uniprot_alignment_pos1, etc}, {label_seq_id2: uniprot_alignment_pos2, etc}]
		self.label_uniprot_dict_pickle_file = osp.join(data_dir, self.name + "_label_uniprot_dict.pickle")
		self.uniprot_label_dict = {} # dict of uniprot -> label mappings of interfacial residues that are mapped between PDB template and UniProt protein, key = ((uniprot1, pdb_template1), (uniprot2, pdb_template2)), value = [{uniprot_alignment_pos1: label_seq_id1, etc}, {uniprot_alignment_pos2: label_seq_id2, etc}]
		self.uniprot_label_dict_pickle_file = osp.join(data_dir, self.name + "_uniprot_label_dict.pickle")

	# get blast alignment
	def get_blast_results(self):
		print('-----Reading BLAST results into a dictionary-----')
		with open(self.blast_file, 'r') as f:
			for line in f:
				items = line.strip().split('\t')
				if items[0] in self.uniprot_pdb_chains:
					if items[1] in self.uniprot_pdb_chains[items[0]]: # only keep the blast results that are in interactions_file
						self.blast_dict.setdefault((items[0], items[1]), []).append(items[2:])

	# get best blast alignment when there is more than one alignment for a given uniprot, pdb chain pair
	# pick the one with the smallest e-value (seems like one with the smallest e-value is always the first entry)
	def get_best_blast_alignment(self):
		print('-----Selecting the best blast alignment for each uniprot, pdb chain pair based on e-value-----')
		for pair in self.blast_dict:
			alignments = self.blast_dict[pair]
			if len(alignments) > 1:
				evalue_list = []
				for alignment in alignments:
					evalue = float(alignment[8])
					evalue_list.append(evalue)
				min_evalue = min(evalue_list)
				min_index = evalue_list.index(min_evalue)
				# print(evalue_list, min_evalue, min_index)
				if pair in self.best_blast_alignment:
					print('already exists!', pair)
				self.best_blast_alignment[pair] = alignments[min_index]
			else:
				if len(alignments) != 1: # length should be equal to 1
					print('ERROR reading in BLAST alignments...')
				else:
					if pair in self.best_blast_alignment:
						print('already exists!', pair)
					self.best_blast_alignment[pair] = alignments[0]
			# print(pair, self.best_blast_alignment[pair])

	def get_residues_mapped(self, uniprot, pdb, chain, alignment, residues_list_auth):
		# query (q) = uniprot protein, subject (s) = pdb chain
		qseq, sseq = alignment[10], alignment[11]
		# qnums and snums are both in label_seq_id format
		qnums, snums = [int(alignment[4]), int(alignment[5])], [int(alignment[6]), int(alignment[7])]
		# ***label_seq_id starts numbering from 1
		residues_list_label = [int(self.auth_label_dict[pdb][chain][auth_seq_id]) for auth_seq_id in residues_list_auth] # convert label residue to int format
		num_residues = len(residues_list_label)
		num_residues_mapped = 0
		residues_mapped = [] # in uniprot numbering system (may or may not correspond to auth_seq_id)
		qindex, sindex = qnums[0]-1, snums[0]-1 # qstart and sstart

		pdb_chain = '_'.join([pdb, chain])

		label_uniprot_dict = {}
		uniprot_label_dict = {}

		# start mapping
		for i in range(len(qseq)):
			if qseq[i] != '-' and sseq[i] != '-': # when 2 amino acids are aligned (substitutions are also accepted)
				qindex += 1
				sindex += 1
				if sindex in residues_list_label:
					num_residues_mapped += 1
					residues_mapped.append(qindex) # save residue num based on uniprot protein
					label_uniprot_dict[str(sindex)] = str(qindex)
					uniprot_label_dict[str(qindex)] = str(sindex)
			elif qseq[i] == '-' and sseq[i] != '-': # when an amino acid on the subject is aligned to a gap on the query
				sindex += 1
			elif qseq[i] != '-' and sseq[i] == '-': # when an amino acid on the query is aligned to a gap on the subject
				qindex += 1
			else: # both are '-'...is impossible
				print('strange alignment...')

		if qindex != qnums[1] and sindex != snums[1]: 
			print('PROBLEM! Somehow indices got messed up...', qseq, sseq, qindex, qnums[1], sindex, snums[1], snums[0])
		else:
			return num_residues_mapped, num_residues, residues_mapped, label_uniprot_dict, uniprot_label_dict

	# returns T/F (whether >= 50% of interfacial residues of one pdb chain map to its aligned uniprot protein)
	def get_mapping_percentage(self, alignment, pair, residues_list_auth):
		uniprot, pdb_chain = pair
		pdb, chain = pdb_chain.split('_')
		# query (q) = uniprot protein, subject (s) = pdb chain
		qseq, sseq = alignment[10], alignment[11]

		# interate through each amino acid in the alignment
		if len(qseq) != len(sseq):
			print('ERROR! BLAST alignment between', uniprot, 'and', pdb + '_' + chain, 'are not the same length!!')
		else:
			num_residues_mapped, num_residues, residues_mapped, label_uniprot_dict, uniprot_label_dict = self.get_residues_mapped(uniprot, pdb, chain, alignment, residues_list_auth)
			
			if num_residues_mapped / num_residues >= 0.5: 
				return True, residues_mapped, label_uniprot_dict, uniprot_label_dict
			else:
				return False, [], {}, {}


	# check whether >= 50% of interfacial residues map to corresponding uniprot protein for both chains
	# amino acid substitutions count, but mapping to a gap on the uniprot protein doesn't count
	def get_confirmed_interactions(self):
		print('-----Finding interactions with >= 50 percent of interfacial residues that map to corresponding uniprot proteins-----')
		uniprot_pairs = []
		with open(self.interactions_file, 'r') as f, open(self.final_interactions_file, 'w') as f_final:
			# next(f)
			header = f.readline().strip().split('\t')
			new_header = header[:-2] + ['Protein1_Interfacial_Residues_Uniprot', 'Protein2_Interfacial_Residues_Uniprot']
			# write header to f_final
			f_final.write('\t'.join(new_header) + '\n')
			for line in f:
				items = line.strip().split('\t')
				# e.g. ('P23511', '4awl_A') and ('Q13952', '4awl_C')
				pair1, pair2 = ('', ''), ('', '')
				if self.name == 'hiunion':
					pair1 = (items[2], items[4] + '_' + items[5])
					pair2 = (items[3], items[4] + '_' + items[6])
				else:
					pair1 = (items[4], items[6] + '_' + items[7])
					pair2 = (items[5], items[6] + '_' + items[8])
				# print(pair1, pair2)
				if pair1 not in self.blast_dict and pair2 not in self.blast_dict:
					print('EROR! Uniprot pairs:', (pair1, pair2), 'is not in the BLAST results file!')
				elif pair1 not in self.blast_dict:
					print('EROR! Uniprot pairs:', pair1, 'is not in the BLAST results file!')
				elif pair2 not in self.blast_dict:
					print('EROR! Uniprot pairs:', pair2, 'is not in the BLAST results file!')
				else:
					# find % of interfacial residues that map to corresponding uniprot protein
					# need to do for both chains
					alignment1 = self.best_blast_alignment[pair1]
					alignment2 = self.best_blast_alignment[pair2]
					residues1_list, residues2_list = [], [] # these residues are in auth format
					if self.name == 'hiunion':
						residues1_list, residues2_list = items[7].split(','), items[8].split(',')
					else:
						residues1_list, residues2_list = items[9].split(','), items[10].split(',')

					# map both pair1 and pair2
					mapped1, uniprot_mapped_residues1, label_uniprot_dict1, uniprot_label_dict1 = self.get_mapping_percentage(alignment1, pair1, residues1_list)
					mapped2, uniprot_mapped_residues2, label_uniprot_dict2, uniprot_label_dict2 = self.get_mapping_percentage(alignment2, pair2, residues2_list)

					# convert residues back to str format
					uniprot_mapped_residues1_str = [str(res) for res in uniprot_mapped_residues1]
					uniprot_mapped_residues2_str = [str(res) for res in uniprot_mapped_residues2]

					# keep interaction only if both residue lists map at least 50%
					if mapped1 and mapped2:
						# update self.label_uniprot_dict and self.uniprot_label_dict
						self.label_uniprot_dict[(pair1, pair2)] = [label_uniprot_dict1, label_uniprot_dict2]
						self.uniprot_label_dict[(pair1, pair2)] = [uniprot_label_dict1, uniprot_label_dict2]
						# write to self.final_interactions_file
						to_write = []
						residues_lists = [','.join(uniprot_mapped_residues1_str)] + [','.join(uniprot_mapped_residues2_str)]
						if self.name == 'hiunion':
							to_write = items[:7] + residues_lists
						else:
							to_write = items[:9] + residues_lists
						f_final.write('\t'.join(to_write) + '\n')
						if (items[2], items[3]) not in uniprot_pairs and (items[3], items[2]) not in uniprot_pairs: # don't need the reverse interaction
							uniprot_pairs.append((items[2], items[3]))
					# else:
					# 	print('Interactions that do not have at least 50 percent of interfacial residues mapped:', pair1, pair2)
		print('Number of PPIs that have at least 50 percent of interfacial residues that map:', len(uniprot_pairs))
		# check for repeated interaction pairs
		for (p1, p2) in uniprot_pairs:
			if (p2, p1) in uniprot_pairs:
				print('repeated!', (p1, p2))
		# dump self.label_uniprot_dict and self.uniprot_label_dict
		pickle_dump(self.label_uniprot_dict, self.label_uniprot_dict_pickle_file)
		pickle_dump(self.uniprot_label_dict, self.uniprot_label_dict_pickle_file)

	def save_best_blast_alignments_to_file(self):
		print('-----Writing best BLAST alignments to file (those that are in the final interactions file)-----')
		pickle_dump(self.best_blast_alignment, self.best_blast_alignment_pickle_file)
		with open(self.best_blast_alignment_file, 'w') as f:
			for pair in self.best_blast_alignment:
				uniprot, pdb_chain = pair
				f.write('\t'.join([uniprot, pdb_chain] + self.best_blast_alignment[pair]) + '\n')


	def build_structural_interactome(self):
		self.get_blast_results()
		self.get_best_blast_alignment()
		self.get_confirmed_interactions()
		self.save_best_blast_alignments_to_file() 

def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed', 'interactome')

	print('-----HI-Union-----')
	hiunion = StructuralInteractome('hiunion', data_dir)
	hiunion.build_structural_interactome()

	print('-----IntAct-----')
	intact = StructuralInteractome('intact', data_dir)
	intact.build_structural_interactome()


if __name__=='__main__':
	main()
