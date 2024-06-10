'''
Categorize missense mutations based on their edgotypes (edgetic, quasi-null, or quasi-wildtype)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import os.path as osp
import os
import math
import numpy as np
from simple_tools import get_missense_info, pickle_dump

class Nums:
	def __init__(self, data_dir, names):
		self.data_dir = data_dir
		self.names = names

	def get_edgotype_quasi_null_wildtype_nums_all(self):
		for name in self.names:
			# print('Finding number of edgetic, quasi-wildtype, and quasi-null mutations for:', name)
			header_dict, missense_info = get_missense_info(osp.join(self.data_dir, name + '_mutation_edgotypes_quasi_null_wildtype.tsv'))
			self.get_edgotype_quasi_null_wildtype_nums(name, missense_info, header_dict)

	def get_avg_blosum_score(self, blosum_scores_list, num):
		if num == 0 and blosum_scores_list == []:
			return 0
		else:
			return sum(blosum_scores_list)/num

	def get_blosum_score_error_bars(self, blosum_scores_list, sqrt_total_num):
		if blosum_scores_list == []:
			return 0
		else:
			return np.std(blosum_scores_list)/sqrt_total_num

	def get_edgotype_quasi_null_wildtype_nums(self, name, missense_info, header_dict):
		print('-----' + name + '-----')
		num_edgetic, num_quasi_wildtype, num_quasi_null = 0, 0, 0
		total_num_mutations = len(missense_info) - 1
		edgetic_blosum_scores, quasi_null_blosum_scores, quasi_wildtype_blosum_scores = [], [], []
		edgotypes_blosum_scores_pickle_file = osp.join(self.data_dir, name + '_edgotypes_blosum_scores.pickle')
		for i in range(1, len(missense_info)):
			# print(header_dict)
			edgotype = missense_info[i][header_dict['edgotype']]
			blosum_score = int(missense_info[i][header_dict['blosum62_score']])
			if edgotype == 'non-edgetic':
				quasi_null_wildtype = missense_info[i][header_dict['quasi_null_wildtype']]
				if quasi_null_wildtype == 'quasi-wildtype':
					num_quasi_wildtype += 1
					quasi_wildtype_blosum_scores.append(blosum_score)
				else:
					if quasi_null_wildtype != 'quasi-null': # sanity check
						print('Something went wrong with:', missense_info[i])
					else:
						num_quasi_null += 1
						quasi_null_blosum_scores.append(blosum_score)
			else: # edgetic mutation
				num_edgetic += 1
				edgetic_blosum_scores.append(blosum_score)
		# avg_blosum_scores = [sum(quasi_wildtype_blosum_scores)/num_quasi_wildtype, sum(edgetic_blosum_scores)/num_edgetic, sum(quasi_null_blosum_scores)/num_quasi_null]
		avg_blosum_scores = [self.get_avg_blosum_score(quasi_wildtype_blosum_scores, num_quasi_wildtype), self.get_avg_blosum_score(edgetic_blosum_scores, num_edgetic), self.get_avg_blosum_score(quasi_null_blosum_scores, num_quasi_null)]
		sqrt_total_num = math.sqrt(total_num_mutations)
		# blosum_scores_error_bars = [np.std(quasi_wildtype_blosum_scores)/sqrt_total_num, np.std(edgetic_blosum_scores)/sqrt_total_num, np.std(quasi_null_blosum_scores)/sqrt_total_num]
		blosum_scores_error_bars = [self.get_blosum_score_error_bars(quasi_wildtype_blosum_scores, sqrt_total_num), self.get_blosum_score_error_bars(edgetic_blosum_scores, sqrt_total_num), self.get_blosum_score_error_bars(quasi_null_blosum_scores, sqrt_total_num)]

		nums_list = [total_num_mutations, num_quasi_wildtype, num_edgetic, num_quasi_null]

		pickle_dump([nums_list, avg_blosum_scores, blosum_scores_error_bars], edgotypes_blosum_scores_pickle_file)
		print('Total number of mutations:', total_num_mutations)
		print('Number of quasi-wildtype mutations:', num_quasi_wildtype, "{:.3f}".format(num_quasi_wildtype/total_num_mutations))
		print('Number of edgetic mutations:', num_edgetic, "{:.3f}".format(num_edgetic/total_num_mutations))
		print('Number of quasi-null mutations:', num_quasi_null, "{:.3f}".format(num_quasi_null/total_num_mutations))

def main():
	data_dir = osp.join(osp.dirname(__file__), '..', 'data', 'processed', 'edgotypes')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	dbsnp_pathogenicities = ['', 'pathogenic']
	names = []

	for interactome in interactomes:
		for data in mutation_data:
			if data == 'dbsnp':
				for dbsnp_pathogenicity in dbsnp_pathogenicities:
					if dbsnp_pathogenicity == '':
						names.append('_'.join([interactome, data, 'missense', 'mutations']))
					else:
						names.append('_'.join([interactome, data, 'missense', 'mutations', 'pathogenic']))
			else: # data == 'clinvar'
				names.append('_'.join([interactome, data, 'missense', 'mutations']))


	# for name in names:
	# 	print(name)

	n = Nums(data_dir, names)
	n.get_edgotype_quasi_null_wildtype_nums_all()

if __name__=='__main__':
	main()
