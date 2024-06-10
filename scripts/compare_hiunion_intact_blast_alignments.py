'''
Combine BLASTP alignments from HI-Union and IntAct interactomes
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
from simple_tools import pickle_load, pickle_dump

# compare selected blast alignments of uniprot & pdb chain pair for interactions in HI-Union & IntAct (since there are overlaps)
# the selected blast alignment for each uniprot & pdb chain pair should be the same between HI-Union & IntAct (if pick blast)
# need for them to be comparable in order to combine modeller alignments

class CombineBlastAlignments:
	def __init__(self, data_dir):
		self.hiunion_blast_alignments = pickle_load(osp.join(data_dir, 'interactome', 'hiunion' + '_blast_best_alignments.pickle')) # key = (protein, pdb_chain), value = [pident, length, midmatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, protein alignment, template_chain alignment]
		self.intact_blast_alignments = pickle_load(osp.join(data_dir, 'interactome', 'intact' + '_blast_best_alignments_physical.pickle')) # key = (protein, pdb_chain), value = [pident, length, midmatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, protein alignment, template_chain alignment]
		self.all_blast_alignments_pickle_file = osp.join(data_dir, 'interactome', 'all_blast_best_alignments.pickle') # key = (protein, pdb_chain), value = [qstart, qend, sstart, send, protein alignment, template_chain alignment]
		self.all_blast_alignments_reduced_info_pickle_file = osp.join(data_dir, 'edgotypes', 'all_blast_best_alignments_reduced_info.pickle') # key = (protein, pdb_chain), value = [qstart, qend, sstart, send, protein alignment, template_chain alignment]
		self.all_blast_alignments = {}
		self.all_blast_alignments_reduced_info = {} # key = (protein, pdb_chain), value = [qstart, qend, sstart, send, protein alignment, template_chain alignment]

	def get_unique_blast_alignments(self):
		overlapping_keys = list(set(self.hiunion_blast_alignments.keys()) & set(self.intact_blast_alignments.keys()))
		hiunion_only_keys = list(set(self.hiunion_blast_alignments.keys()) - set(self.intact_blast_alignments.keys()))
		intact_only_keys = list(set(self.intact_blast_alignments.keys()) - set(self.hiunion_blast_alignments.keys()))
		print('Number of overlapping blast alignments:', len(overlapping_keys))
		print('Number of HI-Union-specific blast alignments:', len(hiunion_only_keys))
		print('Number of IntAct-specific blast alignments:', len(intact_only_keys))

		# getting overlapping alignments
		for key in overlapping_keys:
			if self.hiunion_blast_alignments[key] == self.intact_blast_alignments[key]:
				# print(self.hiunion_blast_alignments[key], self.intact_blast_alignments[key])
				self.all_blast_alignments[key] = self.hiunion_blast_alignments[key]
			else:
				print('Discrepancies between keys!', self.hiunion_blast_alignments[key], self.self.intact_blast_alignments[key])

		# getting HI-Union-specific blast alignments
		for key in hiunion_only_keys:
			self.all_blast_alignments[key] = self.hiunion_blast_alignments[key]

		# getting IntAct-specific blast alignments
		for key in intact_only_keys:
			self.all_blast_alignments[key] = self.intact_blast_alignments[key]

		# pickle self.all_blast_alignments
		pickle_dump(self.all_blast_alignments, self.all_blast_alignments_pickle_file)

		print('Total number of BLAST alignments:', len(self.all_blast_alignments))

	def get_blast_alignment_reduced_info(self):
		for pair in self.all_blast_alignments:
			indices = self.all_blast_alignments[pair][4:8]
			alignments = self.all_blast_alignments[pair][10:]
			self.all_blast_alignments_reduced_info[pair] = indices + alignments
			# print(self.all_blast_alignments[pair], self.all_blast_alignments_reduced_info[pair])
		pickle_dump(self.all_blast_alignments_reduced_info, self.all_blast_alignments_reduced_info_pickle_file)

def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed')

	blast_alignments = CombineBlastAlignments(data_dir)
	blast_alignments.get_unique_blast_alignments()
	blast_alignments.get_blast_alignment_reduced_info()


if __name__=='__main__':
	main()