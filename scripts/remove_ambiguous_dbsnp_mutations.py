'''
Removes ambiguous dbSNP mutations in which the reference and alternative alleles are ambguious in their HGSVC and HGVSG annotations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from simple_tools import get_missense_info, write_missense_info_to_file, pickle_load
import os.path as osp
import os

'''
* Checks ref and alt alleles (from HGVSG) with those in HGVSC (need to do since missense mutations correspond to HGVSC,
  but alt allele frequency corresponds to HGVSG)
* remove mutations in which ref & alt alleles are ambiguous
* also adds the frame of the mutation on the mRNA codons

* if have more than one entry for the same mutation (same nucleotide & amino acid change),
  select the one with the long RefSeq protein sequence
* ALSO NEED TO CHECK THAT ALLELES IN HGVSC ARE VALID (removes ambiguous mutations)
* if (ref_allele == a['ref_allele'] and alt_allele == a['alt_allele']) or (self.get_allele_complement(ref_allele) == a['ref_allele'] and self.get_allele_complement(alt_allele) == a['alt_allele']): # make sure that ref and alt alleles in chromsome & cDNA(mRNA) match or is the complement (they should!)
# PARTICULARLY FOR CHROMOSOME 5
'''

class RemoveAmbiguousMutations:
	def __init__(self, data_dir, mutation_files_list):
		self.data_dir = data_dir
		self.mutation_files_list = mutation_files_list

	def get_mrna_info(self, mutation_cds):
		'''
		example of mutation_cds: c.194G>C
		'''
		# print(mutation_cds)
		to_continue_to_split, alt_allele = mutation_cds.split('>') # split by '>'
		mrna_allele_pos, ref_allele = to_continue_to_split[2:-1], to_continue_to_split[-1]
		# print((ref_allele, mrna_allele_pos, alt_allele))
		return (ref_allele, mrna_allele_pos, alt_allele)        

	def get_allele_complement(self, allele):
		'''
		A pairs with T, C pairs with G
		Use for matching ref and alt alleles from HGVSG with HGVSC in missense variant info
		'''
		if allele == 'A':
			return 'T'
		elif allele == 'T':
			return 'A'
		elif allele == 'C':
			return 'G'
		else: # is G
			if allele != 'G':
				# print('Something went wrong! This is not a valid allele:', allele)
				return -1
			else:
				return 'C'


	def remove_ambiguous_dbsnp_mutations(self, mutation_file):
		print('Looking at:', mutation_file, '...')
		selected_mutations = []
		header_dict, mutations = get_missense_info(mutation_file)
		new_column_names = mutations[0][:15] + ['frame'] + mutations[0][15:]
		# print(new_column_names)
		selected_mutations.append(new_column_names)

		# select mutations
		# create new_header_dict (for selecting mutations)
		new_header_dict = {}
		c = 0
		for column_name in new_column_names:
			new_header_dict[column_name] = c
			c += 1
		
		for mutation in mutations[1:]:
			rs = mutation[header_dict['ID']]
			HGVSG = mutation[header_dict['HGVSG']]
			HGVSC = mutation[header_dict['HGVSC']]
			chrom_ref_allele, chrom_alt_allele = mutation[header_dict['ref_allele']], mutation[header_dict['alt_allele']]
			mrna_ref_allele, _, mrna_alt_allele = self.get_mrna_info(HGVSC.split(':')[1])
			mrna_ref_codon, mrna_alt_codon = mutation[header_dict['mrna_ref_codon']], mutation[header_dict['mrna_alt_codon']]

			# first check mRNA codons
			# test for chrom 3
			frame = str([i for i in range(len(mrna_ref_codon)) if mrna_ref_codon[i] != mrna_alt_codon[i]][0] + 1) # Python indexing starts from 0

			# check frame of mutation on mRNA codon
			if len(frame) != 1: # should have one different nucleotide
				print(mrna_ref_codon, mrna_alt_codon, frame)
			else:
				new_mutation_columns = mutation[:15] + [frame] + mutation[15:]

				# first check whether HGVSG and HGVSC correspond (based on their reference and alternative alleles)
				if (chrom_ref_allele == mrna_ref_allele and chrom_alt_allele == mrna_alt_allele):
					# continue
					selected_mutations.append(new_mutation_columns)
				else:
					mrna_ref_allele_complement = self.get_allele_complement(mrna_ref_allele)
					mrna_alt_allele_complement = self.get_allele_complement(mrna_alt_allele)
					if mrna_ref_allele == -1 or mrna_alt_allele == -1:
						print(HGVSG, HGVSC)
						print(new_mutation_columns)
					elif chrom_ref_allele == mrna_ref_allele_complement and chrom_alt_allele == mrna_alt_allele_complement:
						selected_mutations.append(new_mutation_columns)
					# else: # may occur often...
					# 	print('Not the same!', HGVSG, HGVSC, chrom_ref_allele, chrom_alt_allele)

		return selected_mutations

	def remove_ambiguous_dbsnp_mutations_all(self):
		for mutation_file in self.mutation_files_list:
			selected_mutations = self.remove_ambiguous_dbsnp_mutations(mutation_file)
			new_mutation_file = osp.basename(mutation_file)[:-8] + '.tsv'
			# print(new_mutation_file)
			print(new_mutation_file, 'contains:', len(selected_mutations)-1, 'unambiguous mutations.')
			write_missense_info_to_file(selected_mutations, osp.join(self.data_dir, new_mutation_file))


def main():
	script_dir = osp.dirname(__file__)
	data_dir = osp.join(script_dir, '..', 'data', 'processed', 'mutations')
	dbsnp_unselected_data_dir = osp.join(data_dir, 'dbsnp_unselected')
	mutation_types = ['dbsnp_missense_mutations', 'dbsnp_nonsense_mutations', 'dbsnp_nonstop_mutations']
	clinical_significances = ['', 'pathogenic']
	mutation_files_list = []
	for mutation_type in mutation_types:
		for clinical_significance in clinical_significances:
			if clinical_significance == '':
				mutation_files_list.append(osp.join(dbsnp_unselected_data_dir, '_'.join([mutation_type, 'all.tsv'])))
			else:
				mutation_files_list.append(osp.join(dbsnp_unselected_data_dir, '_'.join([mutation_type, clinical_significance, 'all.tsv'])))


	s = RemoveAmbiguousMutations(data_dir, mutation_files_list)
	s.remove_ambiguous_dbsnp_mutations_all()

if __name__=='__main__':
	main()
		