'''
Processes ClinVar mutations to obtain Mendelian disease-causing missense, nonsense (and nonstop) mutations
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''


import os.path as osp
import pickle
import re
from simple_tools import get_missense_info, pickle_load, pickle_dump, write_missense_info_to_file


class ClinVar():
	def  __init__(self, processed_data_dir, mutation_data_dir, variant_summary_file):
		self.processed_data_dir = processed_data_dir
		self.mutation_data_dir = mutation_data_dir
		self.swissprot_id_list = pickle_load(osp.join(self.processed_data_dir, 'swissprot_ids_list.pickle'))
		self.refseq_prot_gene_name_uniprot_dict = pickle_load(osp.join(self.processed_data_dir, 'refseq_prot_gene_name_uniprot_dict.pickle')) # already checked that these refseq_prots have protein sequences & that the UniProt proteins are SwissProt reviewed
		self.refseq_mrna_refseq_prot_dict = pickle_load(osp.join(self.processed_data_dir, 'refseq_mrna_refseq_prot_dict.pickle'))
		self.refseq_prot_seq_dict = pickle_load(osp.join(self.processed_data_dir, 'refseq_prot_seq_dict.pickle')) # still need
		self.header_dict, self.variant_info = get_missense_info(variant_summary_file) # get header_dict & variant_info; will be updated self.header_dict
		self.pathogenic_snps = []
		self.pathogenic_missense_mutations = []
		self.pathogenic_missense_mutations_file = osp.join(self.mutation_data_dir, 'clinvar_missense_mutations.tsv') # have not matched flanking sequences
		self.gene_chrom_uniprot_pos_alt_res = {} # key = (gene, chrom, uniprot, prot_res_pos), value = alt_res
		# extra files
		self.pathogenic_nonstop_mutations = []
		self.pathogenic_nonstop_mutations_file = osp.join(self.mutation_data_dir, 'clinvar_nonstop_mutations.tsv')
		self.pathogenic_nonsense_mutations = []
		self.pathogenic_nonsense_mutations_file = osp.join(self.mutation_data_dir, 'clinvar_nonsense_mutations.tsv')
		self.pathogenic_synonymous_mutations = []
		self.pathogenic_synonymous_mutations_file = osp.join(self.mutation_data_dir, 'clinvar_synonymous_mutations.tsv')

	# replaces ' ' with '_' in column names in header 
	def add_underscore_to_header(self): 
		# fix self.header_dict
		new_header_dict = {}
		for column_name in self.header_dict:
			old_column_name = column_name
			new_header_dict[column_name.replace(" ", "_")] = self.header_dict[old_column_name]
		self.header_dict = new_header_dict

		# fix first item in self.variant_info
		column_names = self.variant_info[0]
		new_column_names = [column_name.replace(" ", "_") for column_name in column_names]
		self.variant_info[0] = new_column_names 


	def process_variant_summary_file(self):
		'''
		reads in variants from variant_summary_file
		only keeps variants with 'Assembly' = 'GRCh38', 'OriginSimple' = 'germline', 'Type' = 'single nucleotide variant'
								'ClinicalSignificance' = 'Pathogenic', and 'ReviewStatus' = ['criteria provided, multiple submitters, no conflicts'
																							'criteria provided, single submitter', 
																							'practice guideline', 'reviewed by expert panel']
		and saves them to self.pathogenic_snps
		see https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/ for more details regarding review statuses

		'''
		print('-----Processing variant_summary.txt to get pathogenic missense mutations-----')
		num_assembly, num_origin, num_type, num_sig, num_review = 0, 0, 0, 0, 0
		review_statuses_to_keep = ['criteria provided, multiple submitters, no conflicts',
									'criteria provided, single submitter', 
									'practice guideline', 'reviewed by expert panel']
		# print(review_statuses_to_keep)
		self.pathogenic_snps.append(self.variant_info[0]) # need to include header
		for i in range(1, len(self.variant_info)):
			# print(self.variant_info[i])
			if len(self.variant_info[i]) == len(self.header_dict):
				if self.variant_info[i][self.header_dict['Assembly']] == 'GRCh38':
					num_assembly += 1
					if self.variant_info[i][self.header_dict['OriginSimple']] == 'germline':
						num_origin += 1
						if self.variant_info[i][self.header_dict['Type']] == 'single nucleotide variant':
							num_type += 1
							if self.variant_info[i][self.header_dict['ClinicalSignificance']] == 'Pathogenic':
								num_sig += 1
								if self.variant_info[i][self.header_dict['ReviewStatus']] in review_statuses_to_keep: 
								# at least one gold star & no conflicting interpretations
									num_review += 1
									self.pathogenic_snps.append(self.variant_info[i])
			else:
				print('PROBLEM! Somehow this line has', len(self.variant_info[i]), 'lines instead of', len(self.header_dict))

		print(f'Number of variants in GRCh38: {num_assembly}')
		print(f'Number of germline variants in GRCh38: {num_origin}')
		print(f'Number of SNP germline variants in GRCh38: {num_type}')
		print(f'Number of pathogenic SNP germline variants in GRCh38: {num_sig}')
		print(f'Number of goldstar reviewed pathogenic SNP germline variants in GRCh38: {num_review}')

	# adapted from Mohamed's code
	# https://github.com/MohamedGhadie/dispensable_ppi_content/blob/457bedf72716681a1074ac05cbca3eb18afb7599/code/mutation_processing_tools.py#L587
	def to_one_letter_aa(self, aa):
		"""Convert an amino acid three-letter code to one-letter code.
		Args:
			aa (str): amino acid three-letter code.
		Returns:
			str if valid, otherwise '-'
		"""
		# Ter = termination
		# Met = methionine (start)
		# the 20 amino acids + termination codon
		one_letter = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K', 'Trp': 'W', 
					 'Thr': 'T', 'Asn': 'N', 'Pro': 'P', 'Phe': 'F', 'Ala': 'A', 'Gly': 'G', 
					 'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R', 'Met': 'M', 'Val': 'V', 
					 'Glu': 'E', 'Tyr': 'Y', 'Ter': '*'} # removed 'Asx': 'B' and 'Glx': 'Z' (b/c already have Asn and Gln = actual amino acids, and didn't find Glx nor Asx in variant_summary.txt)
		
		aa = aa.title()
		if aa in one_letter:
			return one_letter[aa] 
		else:
			return '-'

	# adapted from Mohamed's code
	# https://github.com/MohamedGhadie/dispensable_ppi_content/blob/457bedf72716681a1074ac05cbca3eb18afb7599/code/mutation_processing_tools.py#L587
	def decompose_protein_snp(self, mut):
		"""Decompose a protein SNP into reference residue, position, and alternative residue.
			Example: p.Gln232Ter
		Args:
			mut (str): mutation to be decomposed.
		Returns:
			tuple
		"""
		snp = mut.strip()
		m = re.match(r"p\.\D{3}\d+\D{3}$", snp)
		if m:
			ref_res, prot_res_pos, alt_res = snp[2:5], snp[5:-3], snp[-3:]
			if ref_res == '-' or alt_res == '_':
				print('These amino acids do not exist somehow!', snp, m)
			# convert amino acid 3 letter code into one letter code
			return self.to_one_letter_aa(ref_res), int(prot_res_pos), self.to_one_letter_aa(alt_res)
		else:
			print('Cannot decompose protein SNP!', snp, m)
			return None, None, None


	def decompose_cdna_mut(self, mut):
		"""Decompose a cDNA mutation into reference allele, position, and alternative allele.
			Example: c.3136G>A
		Args:
			mut (str): mutation to be decomposed.
		Returns:
			tuple
		"""
		to_continue_to_split, alt_allele = mut.split('>') # split by '>'
		mrna_allele_pos, ref_allele = to_continue_to_split[2:-1], to_continue_to_split[-1]
		# print((ref_allele, mrna_allele_pos, alt_allele))
		return (ref_allele, mrna_allele_pos, alt_allele) 

	# adapted from Mohamed's code
	# https://github.com/MohamedGhadie/dispensable_ppi_content/blob/457bedf72716681a1074ac05cbca3eb18afb7599/code/mutation_processing_tools.py#L587
	def decompose_snp_name(self, name):        
		"""Decompose a ClinVar mutation name into several attributes.
			Example: NM_017547.3(FOXRED1):c.694C>T (p.Gln232Ter)
		Args:
			names (str): mutation name to be decomposed.
		Returns:
			gene_name, rna_acc (mRNA accession name), cdna_mut (mutation in cdna format), ref_res (protein wildtype amino acid residue), 
			prot_res_pos (mutation position on protein), alt_res (protein mutation amino acid residue)
		"""
		# added in -*\w* to match gene names that have a hyphen
		m = re.match(r"(\w{2}_\d+\.\d+)*(\(\w*-*\w*\))*(\:)*(c\.[\w\>\-]*)*(\s)*(\(p\.\w*\))*", name.strip())
		if m:
			rna_acc, gene_name, colon, cdna_mut, space, pr_mut = m.groups()
			# remove gene_name from the parentheses
			# if gene_name is None, get it from 'GeneSymbol'
			if gene_name:
				gene_name = gene_name.strip('()')

			# only continue if pr_mut is not None
			if pr_mut: 
				if not gene_name:
					print('No gene name!', name, rna_acc, gene_name, colon, cdna_mut, space, pr_mut)
				# prot_res_pos = mutation position on protein sequence
				ref_res, prot_res_pos, alt_res = self.decompose_protein_snp(pr_mut[1:-1])
				# print(cdna_mut)
				ref_allele, mrna_allele_pos, alt_allele = self.decompose_cdna_mut(cdna_mut)
				# cdna = complementary DNA (only contains conding sequences)
				return gene_name, rna_acc, ref_allele, mrna_allele_pos, alt_allele, ref_res, prot_res_pos, alt_res
			else: # if pr_mut is None
				return gene_name, rna_acc, None, None, None, None, None, None # no pr_mut, so even if have cdna_mut, cannot use anyway

	# updates self.pathogenic_snps with items in decomposed snp names
	def decompose_snp_names(self):
		print('-----Decomposing SNP names to get necessary info------')

		# add swissprot reviewed protein name to self.pathogenic_snps

		# get SNP 'Name'
		name_index = self.header_dict['Name']
		# add in gene_name as well b/c original entries (e.g. from GeneSymbol) may contain more than one gene
		to_add = ['gene_name', 'refseq_mrna', 'ref_allele', 'mrna_allele_pos', 'alt_allele', 'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res', 'uniprot_protein']
		self.pathogenic_snps[0].extend(to_add)
		# add to self.header_dict
		for item in to_add:
			self.header_dict[item] = len(self.header_dict)
		# print(self.header_dict)
		for i in range(1, len(self.pathogenic_snps)):
			snp_name = self.pathogenic_snps[i][name_index]
			gene_name, refseq_mrna, ref_allele, mrna_allele_pos, alt_allele, ref_res, prot_res_pos, alt_res = self.decompose_snp_name(snp_name)

			# find UniProt protein and RefSeq protein acession

			# get refseq protein (needs to exist!)
			refseq_protein = None
			if refseq_mrna in self.refseq_mrna_refseq_prot_dict:
				if self.refseq_mrna_refseq_prot_dict[refseq_mrna] in self.refseq_prot_seq_dict: # still need to check that refseq_protein has an available sequence (else don't keep)
					refseq_protein = self.refseq_mrna_refseq_prot_dict[refseq_mrna]

			# find UniProt protein
			uniprot_protein = None		
			# either from 'OtherIDs' column in variant_summary.txt
			all_ids = self.pathogenic_snps[i][self.header_dict['OtherIDs']]
			if 'UniProtKB' in all_ids:
				id_list = all_ids.split(',')
				for id_ in id_list:
					if 'UniProtKB' in id_:
						if id_.split(':')[1].split('#')[0] in self.swissprot_id_list: # need to check that is SwissProt-reviewed
							# print(id_.split(':')[1].split('#')[0])
							uniprot_protein = id_.split(':')[1].split('#')[0]
							break
			else: # or using RefSeq protein acession and gene name
				if (refseq_protein, gene_name) in self.refseq_prot_gene_name_uniprot_dict:
					uniprot_protein = self.refseq_prot_gene_name_uniprot_dict[(refseq_protein, gene_name)]
			
			# add to self.pathogenic_snps
			items = [gene_name, refseq_mrna, ref_allele, mrna_allele_pos, alt_allele, refseq_protein, ref_res, prot_res_pos, alt_res, uniprot_protein]			
			self.pathogenic_snps[i].extend(items)



	# select missense mutations from self.pathogenic_snps and update self.pathogenic_missense
	def select_missense_mutations(self):
		print('-----Selecting missense mutations------')
		# get list of mrna accessions (for process_ref_seq.py)
		refseq_mrna_list = []
		# write to missense_mutations_file		
		num_missense, num_nonstop, num_nonsense, num_synonymous = 0, 0, 0, 0
		num_missing_info = 0
		self.pathogenic_missense_mutations.append(self.pathogenic_snps[0])
		self.pathogenic_nonsense_mutations.append(self.pathogenic_snps[0])
		self.pathogenic_nonstop_mutations.append(self.pathogenic_snps[0])
		self.pathogenic_synonymous_mutations.append(self.pathogenic_snps[0])
		
		for i in range(1, len(self.pathogenic_snps)):
			# check that there is no missing info
			if None not in self.pathogenic_snps[i]:
				curr_snp = self.pathogenic_snps[i]
				ref_res, alt_res = curr_snp[self.header_dict['ref_res']], curr_snp[self.header_dict['alt_res']]
				gene_name = curr_snp[self.header_dict['gene_name']]
				chrom = curr_snp[self.header_dict['Chromosome']]
				prot_res_pos = curr_snp[self.header_dict['prot_res_pos']]
				if ref_res is not None and alt_res is not None: # just in case, although already checked
					if ref_res != alt_res: # distinguish synonymous & nonsynonymous mutations
						# distinguish between nonsense (alt_res == '*'), nonstop (ref_res == '*'), and missense mutations 
						if alt_res == '*': # nonsense mutation
							num_nonsense += 1
							self.pathogenic_nonsense_mutations.append(self.pathogenic_snps[i])

						elif ref_res == '*': # nonstop mutation
							num_nonstop += 1
							self.pathogenic_nonstop_mutations.append(self.pathogenic_snps[i])

						else:
							num_missense += 1
							self.pathogenic_missense_mutations.append(self.pathogenic_snps[i])

					else: # synonymous mutations
						num_synonymous += 1
						self.pathogenic_synonymous_mutations.append(self.pathogenic_snps[i])
				else:
					print('Something went wrong...', self.pathogenic_snps[i])			
			else:
				# print('Missing info:', self.pathogenic_snps[i])
				num_missing_info += 1


		print('Number of SNPs with missing info:', num_missing_info)
		print('Number of nonsense mutations (removed):', num_nonsense)
		print('Number of nonstop mutations (removed):', num_nonstop)
		print('Number of missense mutations:', num_missense)
		print('Number of synonymous mutations:', num_synonymous)


		# write self.pathogenic_missense_mutations to file
		write_missense_info_to_file(self.pathogenic_missense_mutations, self.pathogenic_missense_mutations_file)
		write_missense_info_to_file(self.pathogenic_nonsense_mutations, self.pathogenic_nonsense_mutations_file)
		write_missense_info_to_file(self.pathogenic_nonstop_mutations, self.pathogenic_nonstop_mutations_file)
		write_missense_info_to_file(self.pathogenic_synonymous_mutations, self.pathogenic_synonymous_mutations_file)


	def process_clinvar_pathogenic_missense_mutations(self):
		self.add_underscore_to_header()
		self.process_variant_summary_file()
		self.decompose_snp_names()
		self.select_missense_mutations()


def main():
	script_dir = osp.dirname(__file__)
	variant_summary_file = osp.join(script_dir, '..', 'data', 'original', 'variant_summary.txt')
	processed_data_dir = osp.join(script_dir, '..', 'data', 'processed')
	mutation_data_dir = osp.join(processed_data_dir, 'mutations')

	clinvar = ClinVar(processed_data_dir, mutation_data_dir, variant_summary_file)
	clinvar.process_clinvar_pathogenic_missense_mutations()


if __name__=='__main__':
	main()



