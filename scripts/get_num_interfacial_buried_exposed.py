'''
Categorize missense mutations based on their locations (IR, buried-NIR, or exposed-NIR)
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import pickle
import os.path as osp
import os
import numpy as np
import math
from simple_tools import pickle_load, pickle_dump, get_missense_info, write_missense_info_to_file

class Location:
	def __init__(self, data_dir, names):
		self.data_dir = data_dir
		self.names = names

	def get_num_interfacial_buried_exposed_all(self):
		for name in self.names:
			header_dict, missense_info = get_missense_info(osp.join(self.data_dir, name + '_mutation_edgotypes_quasi_null_wildtype.tsv'))
			self.get_num_interfacial_buried_exposed(name, header_dict, missense_info)

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

	def get_num_interfacial_buried_exposed(self, name, header_dict, missense_info):
		total_num = 0
		num_interfacial, num_buried, num_exposed = 0, 0, 0
		interfacial_blosum_scores, buried_blosum_scores, exposed_blosum_scores = [], [], []
		mutation_location_blosum_scores_pickle_file = osp.join(self.data_dir, name + '_mutation_location_blosum_scores.pickle')
		print('-----' + name + '-----')
		for i in range(1, len(missense_info)):
			if len(header_dict) != len(missense_info[i]):
				print('oh no!!')
			ir_mutation = missense_info[i][header_dict['ir_mutation']]
			blosum_score = int(missense_info[i][header_dict['blosum62_score']])
			if ir_mutation == '1':
				num_interfacial += 1
				interfacial_blosum_scores.append(blosum_score)
			else: # non-edgetic
				rsa = float(missense_info[i][header_dict['relative_solvent_accessibility']])
				# print(missense_info[i])
				if rsa > 0.25:
					num_exposed += 1
					exposed_blosum_scores.append(blosum_score)
				else:
					num_buried += 1
					buried_blosum_scores.append(blosum_score)
			total_num += 1
		# avg_blosum_scores = [sum(exposed_blosum_scores)/num_exposed, sum(interfacial_blosum_scores)/num_interfacial, sum(buried_blosum_scores)/num_buried]
		avg_blosum_scores = [self.get_avg_blosum_score(exposed_blosum_scores, num_exposed), self.get_avg_blosum_score(interfacial_blosum_scores, num_interfacial), self.get_avg_blosum_score(buried_blosum_scores, num_buried)]
		sqrt_total_num = math.sqrt(total_num)
		# blosum_scores_error_bars = [np.std(exposed_blosum_scores)/sqrt_total_num, np.std(interfacial_blosum_scores)/sqrt_total_num, np.std(buried_blosum_scores)/sqrt_total_num]
		blosum_scores_error_bars = [self.get_blosum_score_error_bars(exposed_blosum_scores, sqrt_total_num), self.get_blosum_score_error_bars(interfacial_blosum_scores, sqrt_total_num), self.get_blosum_score_error_bars(buried_blosum_scores, sqrt_total_num)]
		
		nums_list = [total_num, num_exposed, num_interfacial, num_buried]

		pickle_dump([nums_list, avg_blosum_scores, blosum_scores_error_bars], mutation_location_blosum_scores_pickle_file)
		
		print('Total number of mutations:', total_num)
		print('Number of exposed-NIR mutations:', num_exposed, "{:.2f}".format(num_exposed/total_num))
		print('Number of IR mutations:', num_interfacial, "{:.2f}".format(num_interfacial/total_num))
		print('Number of buried-NIR mutations:', num_buried, "{:.2f}".format(num_buried/total_num))

		

def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed', 'edgotypes')

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

	l = Location(data_dir, names)
	l.get_num_interfacial_buried_exposed_all()


if __name__=='__main__':
	main()

		

