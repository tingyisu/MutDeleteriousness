'''
Processes .json.bz2 file (from https://ftp.ncbi.nih.gov/snp/latest_release/JSON/) for one chromosome one at at time
Writes all non-pathogenic dbSNP missense mutations that have 1000genomes frequencies for the ALT allele into dbsnp_missense_mutations_chr*.tsv
----------------------------------------------
Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
'''

import argparse
import json
import bz2
import re
from simple_tools import pickle_load
import os.path as osp


class ProcessJson:
	def __init__(self, chromosome, data_dir):
		self.chromosome = chromosome
		self.data_dir = data_dir
		self.rs = {}
		self.placements_with_alleles = []
		self.allele_annotations = []
		self.refseq_prot_gene_name_uniprot_dict = pickle_load(osp.join(self.data_dir, 'refseq_prot_gene_name_uniprot_dict.pickle')) # already checked that these refseq_prots have protein sequences & that the UniProt proteins are SwissProt reviewed
		self.refseq_prot_seq_dict = pickle_load(osp.join(self.data_dir, 'refseq_prot_seq_dict.pickle'))
		self.swissprot_id_list = pickle_load(osp.join(self.data_dir, 'swissprot_ids_list.pickle'))
		self.prot_one_letter_to_three_letter = {'C': 'Cys', 'D': 'Asp', 'S': 'Ser', 'Q': 'Gln', 'K': 'Lys', 'W': 'Trp', 'T': 
												'Thr', 'N': 'Asn', 'P': 'Pro', 'F': 'Phe', 'A': 'Ala', 'G': 'Gly', 'I': 'Ile', 
												'L': 'Leu', 'H': 'His', 'R': 'Arg', 'M': 'Met', 'V': 'Val', 'E': 'Glu', 'Y': 'Tyr', '*': 'Ter'}

	def get_mrna_info(self, mutation_cds):
		'''
		example of mutation_cds: c.194G>C
		'''
		# print(mutation_cds)
		to_continue_to_split, alt_allele = mutation_cds.split('>') # split by '>'
		mrna_allele_pos, ref_allele = to_continue_to_split[2:-1], to_continue_to_split[-1]
		# print((ref_allele, mrna_allele_pos, alt_allele))
		return (ref_allele, mrna_allele_pos, alt_allele)        


	def get_chromosome_info(self, mutation_chr):
		'''
		example of mutation_chr: g.1014228G>A
		'''
		# returns ref_allele and alt_allele (position not needed)
		# print(mutation_chr)
		to_continue_to_split, alt_allele = mutation_chr.split('>') # split by '>'
		ref_allele = to_continue_to_split[-1]
		return (ref_allele, alt_allele)

	def get_placements(self):
		'''
		rs refseq_chromosome alleles
		'''
		# checking with mrnas & proteins may remove some repeated missense variants (i.e. same mutation on chromsome but mapped to diff mRNA & protein transcripts)
		self.rs['alleles'] = []  # holder for one or more variant alleles
		self.rs['mrnas'] = [] # holder for one or more HGVSC, use to find missense variants
		self.rs['proteins_info'] = {} # dict with key = (refseq_protein, ref_res, alt_res), value = (prot_res_pos, HGVSP), use to find missense variants
		for alleleinfo in self.placements_with_alleles:

			# has top level placement (ptlp) and assembly info
			placement_annot = alleleinfo['placement_annot']
			# if placement_annot['seq_type'] == 'refseq_chromosome': # only take if is in RefSeq chromosome 'NC_' format (complete genomic molecules, RefSeq genomic 'NG_' is incomplete)
			if alleleinfo['is_ptlp'] and len(placement_annot['seq_id_traits_by_assembly']) > 0:  # get genomic placement and alleles
				# check the assembly name
				if placement_annot['seq_id_traits_by_assembly'][0]['assembly_name']=='GRCh38.p14':
					for a in alleleinfo['alleles']: # get HGVS instead of SPDI 
						# ref_allele & alt_allele are in terms of NC_ (chromsome sequence) HGVSG
						# want to make sure that both are the same as those in NM_ (mRNA transcript) HGVSC
						# decompose HGVSG
						HGVSG = a['hgvs'] # e.g. NC_000001.11:g.1014228G>A
						# if '=' not in HGVSG and '>' in HGVSG: (CANNOT DO THIS OTHERWISE WILL MESS UP ORDER WITH self.allele_annotations)
						refseq_chromosome, mutation_chr = HGVSG.split(':') # refseq_chromosome = NC_000001.11, mutation_chr = g.1014228G>A
						ref_allele, alt_allele = '', ''
						# '=' in HGVSG or 'del' in HGVSG or 'ins' in HGVSG or 'dup' in HGVSG; only '>' is needed to show that is a substitution
						if '>' not in HGVSG: # ref_allele is the same as alt_allele (e.g. {"allele":{"spdi":{"seq_id":"NC_000001.11","position":1014227,"deleted_sequence":"G","inserted_sequence":"G"}},"hgvs":"NC_000001.11:g.1014228="})
							# get ref_allele and alt_allele from spdi
							ref_allele = a['allele']['spdi']['deleted_sequence']
							alt_allele = a['allele']['spdi']['inserted_sequence']
						else:
							ref_allele, alt_allele = self.get_chromosome_info(mutation_chr)
						self.rs['alleles'].append({'ref_allele': ref_allele, 'alt_allele': alt_allele, 'assembly_name': 'GRCh38.p14', 'refseq_chromosome': refseq_chromosome, 'HGVSG': HGVSG})


	def get_refseq_annot_frequencies_clinical_significance_prot_acc(self):
		'''
		rs refseq info
		'''
		unwanted_clinical_sig = ['conflicting-interpretations-of-pathogenicity', 'risk-factor', 'pathogenic', 'likely-pathogenic', 'drug-response']

		# rs_obj['primary_snapshot_data']['placements_with_allele'] and
		# self.allele_annotations = rs_obj['primary_snapshot_data']['allele_annotations']
		# have the same ordering
		for idx in range(0, len(self.rs['alleles'])):
			allele_annotation = self.allele_annotations[idx]['assembly_annotation'][0]
			frequency_entries = self.allele_annotations[idx]['frequency']

			# get frequencies
			for entry in frequency_entries:
				if entry['study_name'] == '1000Genomes' and entry['study_version'] == 2 and allele_annotation['seq_id'] == entry['observation']['seq_id']: # study version 2 = 1000Genomes_30x, study version 1 is the non 30x coverage
					if entry['observation']['deleted_sequence'] == self.rs['alleles'][idx]['ref_allele'] and entry['observation']['inserted_sequence'] == self.rs['alleles'][idx]['alt_allele']:
						alt_freq = entry['allele_count'] / entry['total_count']
						self.rs['alleles'][idx]['alt_frequency'] = alt_freq
						# print(entry)
						# print('---')

			# get clinical significance
			clinical_items = self.allele_annotations[idx]['clinical']
			if clinical_items != []:
				all_clinical_significances = []
				for annot in clinical_items:
					all_clinical_significances.extend(annot['clinical_significances'])

					# get protein accession
					# print(clinical_items)
					for entry in annot['variant_identifiers']:
						organization_info = entry['organization']
						if organization_info == 'UniProtKB':
							uniprot_protein = entry['accession'].split('#')[0]
							if 'uniprot_protein' in self.rs['alleles'][idx]:
								prev_uniprot_protein = self.rs['alleles'][idx]['uniprot_protein']
								if prev_uniprot_protein != uniprot_protein: # if have ambiguity, then find UniProt protein myself using RefSeq protein & gene name
									# print('Something went wrong! Have different UniProtKB protein IDs!', self.rs['alleles'][idx])
									self.rs['alleles'][idx]['uniprot_protein'] = '?'
							else:
								self.rs['alleles'][idx]['uniprot_protein'] = uniprot_protein

				
				# save if don't have unwanted_clinicial_sig
					# print(annot)
				likely_benign = True
				for sig in all_clinical_significances:
					if sig in unwanted_clinical_sig:
						likely_benign = False
				self.rs['alleles'][idx]['likely_benign'] = likely_benign
				self.rs['alleles'][idx]['clinical_significances'] = ','.join(all_clinical_significances)
				if 'uniprot_protein' not in self.rs['alleles'][idx]: # have only one clinical significance annotation & it doesn't come w/ a uniprot protein
					self.rs['alleles'][idx]['uniprot_protein'] = '?'


			else: # have no clinical significance or UniProtKB info
				self.rs['alleles'][idx]['likely_benign'] = True
				self.rs['alleles'][idx]['clinical_significances'] = '?'
				self.rs['alleles'][idx]['uniprot_protein'] = '?'


			# get only RefSeq annotation on NC
			# if self.rs['alleles'][idx]['refseq_chromosome'] == allele_annotation['seq_id']: # removes SOME SNPs (which may be repeated across chromosomes...)
			'''
			GRCh38.p14 chr X    NC_000023.11:g.2738242A>G
			GRCh38.p14 chr X    NC_000023.11:g.2738242A>T
			GRCh38.p14 chr Y    NC_000024.10:g.2738242A>G
			GRCh38.p14 chr Y    NC_000024.10:g.2738242A>T
			'''
			if (re.match('^NC_', allele_annotation['seq_id'])):
				if allele_annotation['genes'] != []:
					self.rs['alleles'][idx]['refseq_annot'] = allele_annotation['genes']

					# print(self.rs['alleles'][idx], len(self.rs['alleles'][idx]['refseq_annot']))



	def process_dbsnp_json_bz2(self):

		# with open('../data/original/test_dbsnp_rs4525.json', 'r') as f:
		#     rs_obj = json.load(f)

		# i = 0
		print('Processing .json.bz2 file for chromosome', self.chromosome, '...')

		# ADD BACK LATER
		# 'refsnp-chr' + self.chromosome + '.json.bz2'
		# 'test_dbsnp', 'rs396190.json.bz2'
		with bz2.BZ2File(osp.join('refsnp-chr' + self.chromosome + '.json.bz2'), 'rb') as f_in, \
			open(osp.join(self.data_dir, 'dbsnp_missense_mutations_chr' + self.chromosome + '.tsv'), 'w') as f_write, \
			open(osp.join(self.data_dir, 'dbsnp_missense_mutations_pathogenic_chr' + self.chromosome + '.tsv'), 'w') as f_write_pathogenic, \
			open(osp.join(self.data_dir, 'dbsnp_nonstop_mutations_chr' + self.chromosome + '.tsv'), 'w') as f_write_nonstop, \
			open(osp.join(self.data_dir, 'dbsnp_nonstop_mutations_pathogenic_chr' + self.chromosome + '.tsv'), 'w') as f_write_nonstop_pathogenic, \
			open(osp.join(self.data_dir, 'dbsnp_nonsense_mutations_chr' + self.chromosome + '.tsv'), 'w') as f_write_nonsense, \
			open(osp.join(self.data_dir, 'dbsnp_nonsense_mutations_pathogenic_chr' + self.chromosome + '.tsv'), 'w') as f_write_nonsense_pathogenic:

			# header line
			f_write.write('\t'.join(['ID', 'ref_allele', 'alt_allele', 'alt_frequency', 'gene_symbol', 'assembly_name', 'chromosome', 
									'HGVSG', 'HGVSC', 'HGVSP', 'refseq_chromosome', 'refseq_mrna', 'mrna_ref_codon', 'mrna_allele_pos', 'mrna_alt_codon', 
									'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res','clinical_significances', 'uniprot_protein']) + '\n')

			f_write_pathogenic.write('\t'.join(['ID', 'ref_allele', 'alt_allele', 'alt_frequency', 'gene_symbol', 'assembly_name', 'chromosome', 
									'HGVSG', 'HGVSC', 'HGVSP', 'refseq_chromosome', 'refseq_mrna', 'mrna_ref_codon', 'mrna_allele_pos', 'mrna_alt_codon', 
									'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res','clinical_significances', 'uniprot_protein']) + '\n')

			f_write_nonstop.write('\t'.join(['ID', 'ref_allele', 'alt_allele', 'alt_frequency', 'gene_symbol', 'assembly_name', 'chromosome', 
									'HGVSG', 'HGVSC', 'HGVSP', 'refseq_chromosome', 'refseq_mrna', 'mrna_ref_codon', 'mrna_allele_pos', 'mrna_alt_codon', 
									'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res','clinical_significances', 'uniprot_protein']) + '\n')

			f_write_nonstop_pathogenic.write('\t'.join(['ID', 'ref_allele', 'alt_allele', 'alt_frequency', 'gene_symbol', 'assembly_name', 'chromosome', 
									'HGVSG', 'HGVSC', 'HGVSP', 'refseq_chromosome', 'refseq_mrna', 'mrna_ref_codon', 'mrna_allele_pos', 'mrna_alt_codon', 
									'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res','clinical_significances', 'uniprot_protein']) + '\n')

			f_write_nonsense.write('\t'.join(['ID', 'ref_allele', 'alt_allele', 'alt_frequency', 'gene_symbol', 'assembly_name', 'chromosome', 
									'HGVSG', 'HGVSC', 'HGVSP', 'refseq_chromosome', 'refseq_mrna', 'mrna_ref_codon', 'mrna_allele_pos', 'mrna_alt_codon', 
									'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res','clinical_significances', 'uniprot_protein']) + '\n')

			f_write_nonsense_pathogenic.write('\t'.join(['ID', 'ref_allele', 'alt_allele', 'alt_frequency', 'gene_symbol', 'assembly_name', 'chromosome', 
									'HGVSG', 'HGVSC', 'HGVSP', 'refseq_chromosome', 'refseq_mrna', 'mrna_ref_codon', 'mrna_allele_pos', 'mrna_alt_codon', 
									'refseq_protein', 'ref_res', 'prot_res_pos', 'alt_res','clinical_significances', 'uniprot_protein']) + '\n')


			for line in f_in:
				rs_obj = json.loads(line.decode('utf-8'))
				# i +=1 

				# if i > 40000:
				#     return
	   
				# re-initialize for the next SNP
				self.rs = {}
				self.placements_with_alleles = []
				self.allele_annotations = []
	   
				self.rs['id'] = rs_obj['refsnp_id']
				if 'primary_snapshot_data' in rs_obj:

					self.placements_with_alleles = rs_obj['primary_snapshot_data']['placements_with_allele']
					self.allele_annotations = rs_obj['primary_snapshot_data']['allele_annotations']
					
					# update self.rs
					self.get_placements()
					self.get_refseq_annot_frequencies_clinical_significance_prot_acc()

					
					for a in self.rs['alleles']:

						if 'refseq_annot' in a:
							for refseq_annot in a['refseq_annot']:
								# print(refseq_annot, '\n')
								# print(a)
								rnas = refseq_annot['rnas']
								gene_symbol = refseq_annot['locus']

								for r in rnas:

									if 'codon_aligned_transcript_change' in r:

										if 'hgvs' in r:
										
											# get mrna info
											HGVSC = r['hgvs']
											
											# print(HGVSC)
											if '=' not in HGVSC and '>' in HGVSC: # otherwise could have e.g. XR_007068461.1:n.5898= 
												refseq_mrna, mutation_cds = HGVSC.split(':')
												ref_allele, mrna_allele_pos, alt_allele = self.get_mrna_info(mutation_cds)
												mrna_ref_codon, mrna_alt_codon = r['codon_aligned_transcript_change']['deleted_sequence'], r['codon_aligned_transcript_change']['inserted_sequence']
												
												if 'protein' in r:
													if 'sequence_ontology' in r['protein']:
														if r['protein']['sequence_ontology'] != []:
															if 'alt_frequency' in a:

																'''
																does not have HGVSP (HGVS for protein) so use SPDI
																HGVSP position is SPDI position + 1 (since SPDI counts the first nucleotide as 0)
																"(P)osition: (a non-negative integer) the number of nucleotides 
																to advance on the reference sequence from the boundary before the 
																first nucleotide, which can be thought of as a 0-based interbase 
																coordinate for the variant start" (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7523648/`)
																'''
																protein = r['protein']['variant']['spdi'] # does not have 'hgvs' entry here (so cannot match HGVSP)
																refseq_protein = protein['seq_id']
																ref_res = protein['deleted_sequence']
																alt_res = protein['inserted_sequence']
																prot_res_pos = str(protein['position'] + 1) # SPDI starts counting from 0
																HGVSP = refseq_protein + ':p.' + self.prot_one_letter_to_three_letter[ref_res] + prot_res_pos + self.prot_one_letter_to_three_letter[alt_res]

																# get SwissProt reviewed uniprot protein
																uniprot_protein = ''
																if a['uniprot_protein'] != '?':
																	if uniprot_protein in self.swissprot_id_list and refseq_protein in self.refseq_prot_seq_dict: # make sure that uniprot protein is swissprot reviewed & that refseq_prot has existing protein sequence
																		uniprot_protein = a['uniprot_protein']
																else: # get uniprot protein from mappings
																	if (refseq_protein, gene_symbol) in self.refseq_prot_gene_name_uniprot_dict: # already checked that uniprot protein is swissprot reviewed & that refseq_prot has existing protein sequence
																		uniprot_protein = self.refseq_prot_gene_name_uniprot_dict[(refseq_protein, gene_symbol)]

																# print(refseq_protein, ref_res, alt_res, self.rs['proteins_info'])

																if uniprot_protein != '': # is a SwissProt reviewed uniprot ID, then write to file

																	# missense mutation
																	if r['protein']['sequence_ontology'][0]['name'] == 'missense_variant':

																		if a['likely_benign']:
																			f_write.write("\t".join(['rs' + self.rs['id'], 
																									a['ref_allele'], 
																									a['alt_allele'], 
																									str(a['alt_frequency']),
																									gene_symbol,
																									a['assembly_name'],
																									self.chromosome, 
																									a['HGVSG'],
																									HGVSC,
																									HGVSP,
																									a['refseq_chromosome'],
																									refseq_mrna,
																									mrna_ref_codon,
																									mrna_allele_pos,
																									mrna_alt_codon,                                                                                              
																									refseq_protein,
																									ref_res,
																									prot_res_pos,
																									alt_res,
																									a['clinical_significances'],
																									uniprot_protein]) + '\n')
																			f_write.flush()

																		else:
																			f_write_pathogenic.write("\t".join(['rs' + self.rs['id'], 
																												a['ref_allele'], 
																												a['alt_allele'], 
																												str(a['alt_frequency']),
																												gene_symbol,
																												a['assembly_name'],
																												self.chromosome, 
																												a['HGVSG'],
																												HGVSC,
																												HGVSP,
																												a['refseq_chromosome'],
																												refseq_mrna,
																												mrna_ref_codon,
																												mrna_allele_pos,
																												mrna_alt_codon,                                                                                                
																												refseq_protein,
																												ref_res,
																												prot_res_pos,
																												alt_res,
																												a['clinical_significances'],
																												uniprot_protein]) + '\n')
																			f_write_pathogenic.flush() 

																	# non-stop mutation
																	elif r['protein']['sequence_ontology'][0]['name'] == 'stop_gained':
																		if a['likely_benign']:
																			f_write_nonsense.write("\t".join(['rs' + self.rs['id'], 
																									a['ref_allele'], 
																									a['alt_allele'], 
																									str(a['alt_frequency']),
																									gene_symbol,
																									a['assembly_name'],
																									self.chromosome, 
																									a['HGVSG'],
																									HGVSC,
																									HGVSP,
																									a['refseq_chromosome'],
																									refseq_mrna,
																									mrna_ref_codon,
																									mrna_allele_pos,
																									mrna_alt_codon,                                                                                               
																									refseq_protein,
																									ref_res,
																									prot_res_pos,
																									alt_res,
																									a['clinical_significances'],
																									uniprot_protein]) + '\n')
																			f_write_nonsense.flush()

																		else:
																			f_write_nonsense_pathogenic.write("\t".join(['rs' + self.rs['id'], 
																												a['ref_allele'], 
																												a['alt_allele'], 
																												str(a['alt_frequency']),
																												gene_symbol,
																												a['assembly_name'],
																												self.chromosome, 
																												a['HGVSG'],
																												HGVSC,
																												HGVSP,
																												a['refseq_chromosome'],
																												refseq_mrna,
																												mrna_ref_codon,
																												mrna_allele_pos,
																												mrna_alt_codon,                                                                                                
																												refseq_protein,
																												ref_res,
																												prot_res_pos,
																												alt_res,
																												a['clinical_significances'],
																												uniprot_protein]) + '\n')
																			f_write_nonsense_pathogenic.flush() 

																	elif r['protein']['sequence_ontology'][0]['name'] == 'stop_lost':
																		if a['likely_benign']:
																			f_write_nonstop.write("\t".join(['rs' + self.rs['id'], 
																									a['ref_allele'], 
																									a['alt_allele'], 
																									str(a['alt_frequency']),
																									gene_symbol,
																									a['assembly_name'],
																									self.chromosome, 
																									a['HGVSG'],
																									HGVSC,
																									HGVSP,
																									a['refseq_chromosome'],
																									refseq_mrna,
																									mrna_ref_codon,
																									mrna_allele_pos,
																									mrna_alt_codon,                                                                                                
																									refseq_protein,
																									ref_res,
																									prot_res_pos,
																									alt_res,
																									a['clinical_significances'],
																									uniprot_protein]) + '\n')
																			f_write_nonstop.flush()

																		else:
																			f_write_nonstop_pathogenic.write("\t".join(['rs' + self.rs['id'], 
																												a['ref_allele'], 
																												a['alt_allele'], 
																												str(a['alt_frequency']),
																												gene_symbol,
																												a['assembly_name'],
																												self.chromosome, 
																												a['HGVSG'],
																												HGVSC,
																												HGVSP,
																												a['refseq_chromosome'],
																												refseq_mrna,
																												mrna_ref_codon,
																												mrna_allele_pos,
																												mrna_alt_codon,                                                                                               
																												refseq_protein,
																												ref_res,
																												prot_res_pos,
																												alt_res,
																												a['clinical_significances'],
																												uniprot_protein]) + '\n')
																			f_write_nonstop_pathogenic.flush()

def main():

	parser = argparse.ArgumentParser()
	parser.add_argument('-c', '--chromosome') # '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y'
	parser.add_argument('-d', '--data_dir') # /home/user/scratch/username/dbsnp/data (replace username with your compute canada username)
	args = parser.parse_args()


	p = ProcessJson(args.chromosome, args.data_dir)
	p.process_dbsnp_json_bz2()
	
   

if __name__=='__main__':
	main()