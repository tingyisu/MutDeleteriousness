'''
Removes redundant mutations with the same amino acid change at a given residue position
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

from simple_tools import get_missense_info, write_missense_info_to_file, pickle_load
import os.path as osp
import os

# keeps only one mutation with the same amino acid change (even if diff nucleotide change) at a given position on a UniProt protein 
# and for dbSNP mutations specifically, which includes the same mutation on diff mRNA/protein transcripts, keeps only one (with the longest protein transcript)

class RemoveRedundant:
	def __init__(self, data_dir, final_data_dir, refseq_prot_seq_dict_pickle_file, mutation_files_list):
		self.data_dir = data_dir
		self.final_data_dir = final_data_dir
		self.refseq_prot_seq_dict = pickle_load(refseq_prot_seq_dict_pickle_file)
		if not osp.exists(self.final_data_dir):
			os.mkdir(self.final_data_dir)
		self.mutation_files_list = mutation_files_list

	# for dbSNP mutations only
	def find_longest_protein_transcript(self, proteins_list):
		protein_lengths = [len(self.refseq_prot_seq_dict[protein]) for protein in proteins_list]
		# print(protein_lengths)
		return protein_lengths.index(max(protein_lengths)) # if more than one max, returns the first occurrence

	def select_dbsnp_mutation(self, mutation_file):
		header_dict, mutations = get_missense_info(mutation_file)
		print('# of mapped mutations:', len(mutations)-1)

		# get only one mutation per (rs, chrom_ref_allele, chrom_alt_allele)
		mutations_dict = {} # key = (rs, chrom_ref_allele, chrom_alt_allele), value = list of mutations
		for mutation in mutations[1:]:
			rs = mutation[header_dict['ID']]
			chrom_ref_allele = mutation[header_dict['ref_allele']]
			chrom_alt_allele = mutation[header_dict['alt_allele']]
			if (rs, chrom_ref_allele, chrom_alt_allele) not in mutations_dict:
				mutations_dict[(rs, chrom_ref_allele, chrom_alt_allele)] = [mutation]
			else:
				mutations_dict[(rs, chrom_ref_allele, chrom_alt_allele)].append(mutation)

		# sanity check
		# print(len(mutations_dict.keys()))

		# for the same mutation on diff mRNA/protein transcripts, keeps only one (with the longest protein transcript)
		selected_mutations = [mutations[0]]
		for key in mutations_dict:
			proteins_list = [unambiguous_mutation[header_dict['refseq_protein']] for unambiguous_mutation in mutations_dict[key]]
			selected_mutation = mutations_dict[key][self.find_longest_protein_transcript(proteins_list)]
			selected_mutations.append(selected_mutation)
		return header_dict, selected_mutations

	def remove_redundant_mutations(self, header_dict, mutations):
		# only keep one mutation per (uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res)
		mutations_dict = {} # key = (uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res), value = list of mutations
		for mutation in mutations[1:]:
			uniprot_protein = mutation[header_dict['uniprot_protein']]
			prot_res_pos_in_uniprot_protein = mutation[header_dict['prot_res_pos_in_uniprot_protein']]
			ref_res = mutation[header_dict['ref_res']]
			alt_res = mutation[header_dict['alt_res']]
			if (uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res) not in mutations_dict:
				mutations_dict[(uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res)] = [mutation]
			else:
				mutations_dict[(uniprot_protein, prot_res_pos_in_uniprot_protein, ref_res, alt_res)].append(mutation)

		# sanity check
		# print(len(mutations_dict.keys()))

		# take first instance 
		selected_mutations = [mutations[0]]
		for key in mutations_dict:
			selected_mutations.append(mutations_dict[key][0]) # take first one

		return selected_mutations


	def remove_redundant_mutations_all(self):
		for mutation_file in self.mutation_files_list:
			mutation_file_path = osp.join(self.data_dir, mutation_file)
			nonredundant_mutations = []
			if 'dbsnp' in mutation_file:
				print('-----', mutation_file, '-----')
				
				header_dict, dbsnp_selected_mutations = self.select_dbsnp_mutation(mutation_file_path)
				nonredundant_mutations = self.remove_redundant_mutations(header_dict, dbsnp_selected_mutations)
				print('# of selected mutations (one per rs, chrom_ref_alleles, chrom_alt_allele):', len(dbsnp_selected_mutations)-1)
				print('# of nonredundant mutations:', len(nonredundant_mutations)-1)

				mutation_file_items = mutation_file.split('_')[:-1]
				write_missense_info_to_file(nonredundant_mutations, osp.join(self.final_data_dir, '_'.join(mutation_file_items + ['mutations', 'nonredundant.tsv'])))
			else:
				print('-----', mutation_file, '-----')
				header_dict, mutations = get_missense_info(mutation_file_path)
				nonredundant_mutations = self.remove_redundant_mutations(header_dict, mutations)
				print('# of mapped mutations:', len(mutations)-1)
				print('# of nonredundant mutations:', len(nonredundant_mutations)-1)

				mutation_file_items = mutation_file.split('_')[:-1]
				write_missense_info_to_file(nonredundant_mutations, osp.join(self.final_data_dir, '_'.join(mutation_file_items + ['mutations', 'nonredundant.tsv'])))

def main():
	script_dir = osp.dirname(__file__)
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')
	data_dir = osp.join(processed_data_dir, 'mutations')
	final_data_dir = osp.join(script_dir, '..', 'data', 'processed', 'mutations_final')
	refseq_prot_seq_dict_pickle_file = osp.join(processed_data_dir, 'refseq_prot_seq_dict.pickle')

	interactomes = ['hiunion', 'intact']
	mutation_data = ['clinvar', 'dbsnp']
	mutation_types = ['missense', 'nonstop', 'nonsense']
	mutation_files_list = []

	for interactome in interactomes:
		for data in mutation_data:
			for mutation_type in mutation_types:
				mutation_files_list.append('_'.join([interactome, 'mapped', data, mutation_type, 'mutations.tsv']))
							
	r = RemoveRedundant(data_dir, final_data_dir, refseq_prot_seq_dict_pickle_file, mutation_files_list)
	r.remove_redundant_mutations_all()


if __name__=='__main__':
	main()
