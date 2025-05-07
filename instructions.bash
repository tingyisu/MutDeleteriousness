#!/bin/bash

# '''
# Instructions for running scripts in MutDeleteriousness
# ----------------------------------------------
# Author: Ting-Yi Su (ting-yi.su@mail.mcgill.ca)
# '''

# Please follow the steps below carefully
# Some steps are very memory/time consuming to run and should ideally be run in parallel on a server (scripts for submitting SLURM jobs are provided for these steps).

# first make sure you're in the MutDeleteriousness directory
# and create the following folders to store data
mkdir data
cd data
mkdir original
mkdir processed

# make subfolders in the data/processed dir
cd processed
mkdir interactome
mkdir mutations
mkdir mutations_final
mkdir edgotypes
mkdir edgotypes_final
mkdir nonsense_on_si


# ----------DOWNLOAD THE FOLLOWING DATA FILES:----------
cd ../original
# **1. UniProt mappings file**
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
gunzip HUMAN_9606_idmapping.dat.gz

# **2. RefSeqGene mappings (for mutations; mapping mRNA accession NM_* to protein accession NP_*)** 
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene

# **3. RefSeq protein sequence files**
mkdir ref_seq
cd ref_seq
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.*.protein.faa.gz
gunzip human.*.protein.faa.gz
cat * > ../../processed/merged_ref_seq_protein.fa

# **4. Uniprot reviewed (UniProtKB/SWISS-PROT) human reference proteome in list and FASTA format**
cd ../
wget 'https://rest.uniprot.org/uniprotkb/stream?format=list&query=%28%28taxonomy_id%3A9606%29%20AND%20%28reviewed%3Atrue%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28reviewed%3Atrue%29' -O uniprot_reviewed_human_proteome.list
wget 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%20AND%20%28reviewed%3Atrue%29%29%20AND%20%28model_organism%3A9606%29%20AND%20%28reviewed%3Atrue%29' -O uniprot_reviewed_human_proteome.fasta

# **5. HI-Union human binary protein-protein interaction dataset**
wget http://www.interactome-atlas.org/data/HI-union.tsv

# **6. IntAct human binary protein-protein interaction dataset**
wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip
unzip intact.zip
rm intact.zip intact_negative.txt

# **7. Protein Data Bank (PDB) SeqRes chain sequences**
wget https://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt

# **8. File containing the resolutions of PDB structures**
wget ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx

# **9. ClinVar (Mendelian disease-causing) mutations**
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
gunzip variant_summary.txt.gz


# ----------SELECT PDB STRUCTURAL TEMPLATES FOR HI-UNION AND INTACT REFERENCE PPIS----------
cd ../../scripts
# **1. Process all mappings necessary for processing & mapping mutations to the structural interactomes**
python3 process_mappings.py

# **2. Process HI-union and IntAct data to get human reference PPIs**
python3 process_hiunion_data.py # also processes PDB SeqRes sequences for input into BLASTP
python3 process_intact_data.py

# **3. Perform BLASTP sequence alignment of SwissProt proteins (within HI-union and IntAct) against PDB chains**
cd ../data/processed/interactome
makeblastdb -in pdb_seqres_blast.fasta -dbtype prot -out pdb_seqres_db
# -----HI-union-----
blastp -db pdb_seqres_db -query hiunion_uniprot_sequences_blast.fasta -out hiunion_blast_results_with_seq.tsv -outfmt "6 std qseq sseq"
# keep all alignments with E-value < 10^-5 and pick the best alignment (listed first) for each pair of SwissProt protein and PDB chain
awk -F"\t" '$11<10^-5' "hiunion_blast_results_with_seq.tsv" > "hiunion_blast_results_with_seq_eval-5.tsv"
cut -d$'\t' -f 1-2 hiunion_blast_results_with_seq_eval-5.tsv | awk -F"\t" '!seen[$1, $2]++' > hiunion_blast_for_finding_interactions.tsv
# -----IntAct-----
blastp -db pdb_seqres_db -query intact_uniprot_sequences_blast.fasta -out intact_blast_results_with_seq.tsv -outfmt "6 std qseq sseq"
# keep all alignments with E-value < 10^-5 and pick the best alignment (listed first) for each pair of SwissProt protein and PDB chain
awk -F"\t" '$11<10^-5' "intact_blast_results_with_seq.tsv" > "intact_blast_results_with_seq_eval-5.tsv"
cut -d$'\t' -f 1-2 intact_blast_results_with_seq_eval-5.tsv | awk -F"\t" '!seen[$1, $2]++' > intact_blast_for_finding_interactions.tsv

# **4. Use a distance-based approach to find interfacial residues between interacting chain-pairs in PDB structures**
# NOTE: Please use a server for this step as it will take several months if you don't run multiple jobs in parallel on a server 
# The server we used was Compute Canada/Digital Research Alliance of Canada, which uses the SLURM job scheduler

# (1) CREATE A SCRIPT DIRECTORY called 'interfacial_residues' on your server and COPY the following scripts over (this is where you will run your scripts):
# If you're using Compute Canada, you should create the 'interfacial_residues' folder within your scratch directory (e.g. /home/username/scratch/interfacial_residues)
# The scratch directory can store up to 20TB of data, but they will be deleted after 60 days
# (1) process_blast_results.py
# (2) get_hiunion_pdb_chain_dict.py
# (3) get_intact_pdb_chain_dict.py
# (4) split_pdb_chain_dict.py
# (5) create_residues_slurm_files.py
# (6) create_bash_for_submitting_residues_slurm_jobs.py
# (7) hiunion_residues_split.py
# (8) intact_residues_split.py
# (9) memo_residues.py
# (10) simple_tools.py
# (11) hiunion_write_interactions.py
# (12) intact_write_interactions.py
# (13) get_auth_label_dict.py
# (14) get_auth_label_dict.slurm

# (2) CREATE AN OUTPUT DIRECTORY (-o or --output_directory) called 'files' under your 'interfacial_residues' folder on your server and COPY the following files from data/processed/interactome there:
# If you're using Compute Canada, you should create the 'interfacial_residues/files' folder within your scratch directory (e.g. /home/username/scratch/interfacial_residues/files)
# (1) hiunion_blast_for_finding_interactions.tsv
# (2) HI-union_unique_uniprot_pairs.pickle
# (3) intact_blast_for_finding_interactions.tsv
# (4) intact_uniprot_pairs_physical.pickle
# (5) HI-union_interacting_pairs.pickle
# (6) intact_gene_dict_physical.pickle
# (7) intact_id_dict_physical.pickle

# (3) CREATE A PDB DIRECTORY (-d or --pdb_download_path) called 'pdb_cif' to store downloaded PDB mmCIF files
# If you're using Compute Canada, you should create the 'pdb_cif' folder within your scratch directory (e.g. /home/username/scratch/pdb_cif)

# the following arguments are REQUIRED (please use ABSOLUTE PATHS):
# -o or --output_directory (this is the output dir you created in step (2))
# -d or --pdb_download_path (this is the pdb dir you created in step (3))
# -u or --num_splits_for_hiunion (how many splits/jobs to create for HI-Union)
# -i or --num_splits_for_intact (how many splits/jobs to create for IntAct)
# -n or --num_days (max number of days to run each job)
# -e or --email_address (input your email address if you would like to receive emails on when a job has completed or failed)
# -a or --account (account name/compute canada group, e.g. def-*)

# run the following commands one-by-one in the script directory on your server:

# (4) Find possible homologs from BLASTP alignment results
python3 process_blast_results.py -o <output_directory> 

# (5) Get HI-union and IntAct pdb_chain_dicts (dictionaries of pairs of PDB chains) & download PDB structures in .cif (mmCIF) format
python3 get_hiunion_pdb_chain_dict.py -d <pdb_download_path> -o <output_directory>
python3 get_intact_pdb_chain_dict.py -d <pdb_download_path> -o <output_directory>

# (6) Split the HI-union pdb_chain_dict into 20 jobs and the IntAct pdb_chain_dict into 50 jobs
python3 split_pdb_chain_dict.py -o <output_directory> -u 20 -i 50 # 20 splits for HI-Union, 50 for IntAct

# (7) Create SLURM job files to find interfacial residues between pairs of PDB chains
# 20 splits for HI-Union, 50 for IntAct, run each job for 5 days
python3 create_residues_slurm_files.py -d <pdb_download_path> -o <output_directory> -u 20 -i 50 -n 5 -e <email_address> -a <account> 
python3 create_bash_for_submitting_residues_slurm_jobs.py -u 20 -i 50 

# (8) Submit SLURM jobs
bash submit_residues_slurm_jobs.bash

# (9) If any of the SLURM jobs were not completed in the specified amount of time, 'JOB FAILED' messages will be sent to your email address
# Simply submit the uncompleted jobs again and they will continue to run from where they stopped at; resubmit as many times as needed
sbatch _name_of_uncompleted_job.slurm # replace _name_of_uncompleted_job with the name of the job (e.g. hiunion_1)

# (10) Once all SLURM jobs are completed, combine the memo files in your -o or --output_directory
cd output_directory # replace with your <output_directory>
cat hiunion_memo_residues_split*.tsv > hiunion_memo_residues_combined.tsv
cat intact_memo_residues_split*.tsv > intact_memo_residues_combined.tsv # excluding overlapping HI-union chain pairs
cat hiunion_memo_residues_combined.tsv intact_memo_residues_combined.tsv > memo_residues_combined.tsv # all chain pairs, need to use for IntAct

# (11) Run the following python scripts to write interfacial residues found in *_memo_residues_combined.tsv to file
cd ..
python3 hiunion_write_interactions.py -o <output_directory>
python3 intact_write_interactions.py -o <output_directory>

# (12) Create a dictionary storing PDB auth_seq_id to label_seq_id mappings for converting auth residues to label residues
sbatch get_auth_label_dict.slurm # change the -d (--pdb_download_path), -o (--output_directory), and sbatch --mail_user and --account arguments

# (13) COPY the following files back to the local ../data/processed/interactome directory, and exit out of your server
# 1. hiunion_uniprot_pdb_chains.pickle
# 2. intact_uniprot_pdb_chains.pickle
# 3. hiunion_auth_label_dict.pickle
# 4. intact_auth_label_dict.pickle
# 5. hiunion_interfacial_residues.tsv
# 6. intact_interfacial_residues_physical.tsv
# and the following additional (optional) files for backup...
# 7. hiunion_memo_residues_combined.tsv
# 8. intact_memo_residues_combined.tsv

# **5. Construct the HI-union and IntAct structural interactomes**

# (1) Create final interactions file (keep interactions with >=50% of interfacial residues mapping to their corresponding Uniprot proteins)
python3 build_structural_interactome.py

# (2) Select PDB structural template for each binary interaction
# threshold for templates: 0.0 < PDB resolution <= 3.5 & BLASTP alignment bitscore >= 50
python3 select_pdb_structural_template.py


# ----------PROCESS AND MAP MISSENSE AND NONSENSE MUTATIONS ONTO STRUCTURAL TEMPLATES----------

# **1. Process ClinVar mutations**
python3 process_clinvar_mutations.py

# **2. Download and process dbSNP mutations on your server**
# NOTE: need to run on a server because dbSNP data is very large (each chromosome file is tens of GB)

# (1) Create a 'dbsnp' directory on your server; all downloaded & processed dbSNP data will be stored here
# If you're using Compute Canada, create the directory your scratch directory as the dbSNP data is hundreds of GB total
# replace, 'username' with your compute canada username
cd /home/username/scratch/
mkdir dbsnp
cd dbsnp

# (2) Copy the following script files from your local script directory to your 'dbsnp' directory on your server
# 1. download_dbsnp_json.slurm
# 2. process_dbsnp_json.py
# 3. simple_tools.py 

# (3) Copy all 'dbsnp_json_*.slurm' files under the local script/dbsnp_json dir to your 'dbsnp' dir on your server

# (4) Create a new 'data' dir under the 'dbsnp' dir on your server and copy the following data files from your local data/processed/ dir:
mkdir data
cd data
# 1. swissprot_ids_list.pickle
# 2. refseq_prot_gene_name_uniprot_dict.pickle
# 3. refseq_prot_seq_dict.pickle

# (5) Download dbSNP data in .json format (one .json file for each chromosome) 
# (need to change the directory names in the .slurm job)
cd ../
sbatch download_dbsnp_json.slurm # change the sbatch --mail_user and --account arguments

# (6) Process the dbSNP .json files
# Change the the -d (--data_dir) and sbatch --mail_user and --account arguments in each .slurm job
# These jobs will process the dbSNP .json files (4 chromsomes per job)
for i in {1..6}
do
	sbatch dbsnp_json_"$i".slurm
done

# (7) Combine the dbSNP output files
# missense mutations
cat data/dbsnp_missense_mutations_chr1.tsv > data/dbsnp_missense_mutations_all.tsv # want the header line at the beginning of the file
# nonsense mutations
cat data/dbsnp_nonsense_mutations_chr1.tsv > data/dbsnp_nonsense_mutations_all.tsv # want the header line at the beginning of the file
# nonstop mutations
cat data/dbsnp_nonstop_mutations_chr1.tsv > data/dbsnp_nonstop_mutations_all.tsv # want the header line at the beginning of the file
chroms=('2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' 'X' 'Y')
for i in "${chroms[@]}"
do
	# missense mutations
	awk FNR!=1 data/dbsnp_missense_mutations_chr"$i".tsv >> data/dbsnp_missense_mutations_all.tsv # skip header line
	# nonsense mutations
	awk FNR!=1 data/dbsnp_nonsense_mutations_chr"$i".tsv >> data/dbsnp_nonsense_mutations_all.tsv # skip header line
	# nonstop mutations
	awk FNR!=1 data/dbsnp_nonstop_mutations_chr"$i".tsv >> data/dbsnp_nonstop_mutations_all.tsv # skip header line
done

# (8) Copy the output files back to your local data/processed/mutations directory
# create a 'dbsnp_unselected' folder under data/processed/mutations
cd data/processed/mutations
mkdir dbsnp_unselected
cd dbsnp_unselected
# move the *all.tsv files on compute canada to dbsnp_unselected
# use scp to move; replace username with your compute canada username
scp username@graham.computecanada.ca:/home/username/scratch/dbsnp/data/*all.tsv .

# (9) Remove ambiguous dbSNP mutations
# i.e. mutations where the HGVSC (nucleotide change on mRNA/cDNA) doesn't correspond with the HGVSG (nucleotide change on chromosome)
# mutations are written to data/processed/mutations
python3 remove_ambiguous_dbsnp_mutations.py

# **3. Map dbSNP and ClinVar mutations to UniProt proteins**
# obtain the mutation flanking sequences (first 10 residues & last 10 residues flanking the mutation)
# then match mutations to uniprot sequences via their flanking sequences
# kept mutations that have diff positions on RefSeq and UniProt
python3 map_mutation_flanking_seq_to_uniprot.py

# **4. Remove redundant mutations**
# for ClinVar & dbSNP mutations, keep only one mutation with the same amino acid change (even if diff nucleotide change) at a given position on a UniProt protein 
# for dbSNP mutations only: for the same mutation on different mRNA/protein transcripts, pick one on longest protein transcript, 
# or if they're the same length, pick the first one encountered
python3 remove_redundant_mutations.py

# **5. Map nonsense mutations onto the HI-union and IntAct structural interactomes**
python3 get_nonsense_mutations_on_si.py 

# **6. Find missense mutations in interfacial residues (IR)**
# also removes missense mutations on proteins that do not participate in protein-protein interactions within the interactomes,
# as well as those that do not fall within BLASTP alignments
python3 get_ir_mutations.py

# ----------BUILD HI-UNION AND INTACT STRUCTURAL INTERACTOMES USING MODELLER----------

# **1. Combine BLASTP alignments from HI-Union and IntAct interactomes**
# the overlapping uniprot, pdb chain pairs should have the same BLASTP alignment
# will need this to run MODELLER (also reduces the number of computations needed)
python3 compare_hiunion_intact_blast_alignments.py

# **2. Obtain all structural templates to be modelled**
# get all pairs of uniprot, pdb chains (structural templates) that need to be modelled
# also prepares IR & non-IR mutations for FoldX energy calculations
python3 get_modeller_templates_foldx_mutations.py 

# **3. Download and run MODELLER** 
# NOTE: Please use a server for this step

# (1) Install MODELLER (https://salilab.org/modeller/10.3/release.html#deb, choose either conda or Linux distributions if using compute canada)
# if you have Anaconda, you can install MODELLER as follows:
conda install -c salilab modeller

# (2) CREATE A SCRIPT DIRECTORY called 'modeller' on your server and move the following local scripts to this directory
# If you're using compute canada, create the 'modeller' folder under your projects directory (e.g. /home/username/projects/def-*/username/modeller) 
# 	1. convert_blast_alignment_to_pdb_numbering.py
#	2. split_modeller_alignments.py
#	3. run_modeller_all.py
#	4. list_ali_files_with_discrepancies_to_fix.py
#	5. run_modeller_uncompleted.py
#	6. prepare_modeller_pdb_files_for_foldx.py
#	7. prepare_modeller_pdb_files_for_foldx.slurm
#	8. get_edgetic_mutations.py
# 	9. get_edgetic_mutations.slurm
#	10. run_modeller_additional.py
#	11. get_additional_foldx_buildmodel_mutations.py
#	12. get_additional_foldx_buildmodel_mutations.slurm
#	13. simple_tools.py
#	14. calc_rsa.py
#	15. get_rsa_slurm_file.py
#	16. get_quasi_null_wildtype.py
#	17. get_quasi_null_wildtype.slurm

# (3) CREATE A DATA DIRECTORY called 'files' folder under the 'modeller' directory on your server
# If you're using compute canada, create the 'modeller/files' folder under your projects directory (e.g. /home/username/projects/def-*/username/modeller/files) 
# Copy the following files in /data/processed/interactome onto compute canada under the 'modeller/files' folder:
#	1. all_blast_best_alignments.pickle
#	2. *_selected_pdb_structural_template.pickle (* = hiunion & intact)
# Copy the following files in (/data/processed/edgotypes) onto compute canada under the 'files' folder:
# 	1. all_blast_best_alignments_reduced_info.pickle
#	2. pdb_structures_to_download.pickle
#	3. all_modeller_templates.pickle
#	4. all_foldx_pssm_mutations.pickle
#	5. all_foldx_buildmodel_mutations.pickle
#	6. *_ir_with_modeller_foldx_mutations.tsv

# The command line arguments for the scripts in step (2) are as follows (please use FULL ABSOLUTE paths):
# 	-s or --script_dir (this is the absolute path of your 'modeller' directory where the scripts in step (2) are under)
#	-d or --data_dir (this is the absolute path for the 'modeller/files' folder where the files are stores in step (3))
#	-p or --pdb_download_dir (this the wherOe you will be downloading your .cif files (if you don't have them already); should preferably be in scratch unless you have enough space elsewhere)
#	-e or --email_address (this is your email address that you want SLURM to send the progress of your jobs to)
#	-a or --account (this is your Compute Canada project resource allocation account/group name, starts with def-*)
#	-c or --scratch_dir (this is your scratch directory on Compute Canada; should be /home/username/scratch)

# (4) Convert BLASTP alignments to correspond with existing atoms in the mmCIF files (some atoms are not structurally resolved)
python3 convert_blast_alignment_to_pdb_numbering.py -d <data_dir> -p <pdb_download_dir>

# (5) create and split .ali alignment files up to run multiple MODELLER jobs in parallel
python3 split_modeller_alignments.py -s <script_dir> -p <pdb_download_dir> -e <email_address> -a <account>

# (6) Submit the run_modeller_part_*.slurm jobs created by split_modeller_alignments.py
# replace n with the number of the last job (will be printed when you run split_modeller_alignments.py)
for i in {0..n}
do
	sbatch run_modeller_part_"$i".slurm
done

# (7) Remove unnecessary intermediate files
# run the following to keep only .pdb files and move them to a 'results' directory
DIR=_insert_your_script_directory # e.g DIR=/home/username/projects/def-*/username/modeller
rm "$DIR"/*.ini "$DIR"/*.rsr "$DIR"/*.sch "$DIR"/*.D00000001 "$DIR"/*.V99990001
mkdir results
mv $DIR/*.pdb "$DIR"/results

# (8) Resubmit any uncompleted SLURM jobs
# if any of the run_modeller_part_"$i".slurm were not completed in the allocated time, simply submit the .slurm job again
sbatch run_modeller_part_*.slurm

# (9) Fix the alignment files listed when running the following script:
python3 list_ali_files_with_discrepancies_to_fix.py -s <script_dir>
# These alignments have failed to produce MODELLER .pdb files or produced .pdb files with missing coordinates because:
# 1. index issues in .ali file, possibly due to having a 'CRO' heteroatom in the mmCIF file, where the corresponding residues are 'TYG' (3 residues, so changes the indices by 2) in the PDB SeqRes and consequently the alignment also
#    --> to fix, -2 in the start & end indices for the structure part in the .ali file
# 	 --> occurs particularly for 6wvg_A and 6wvg_B chains
# e.g. >P1;6wvg_B
# structure:/home/username/scratch/pdb_cif/6wvg.cif:153:B:357:B::::
#KLLKYVLFFFNLLFWICGCCILGFGIYLLIHNNFGVLF----HNLPSLTLGNVFVIVGS--IIMVVAFLGCMGSIKENKSLLMSFFILLLIILLAEVTLAILLFVYEQKLNEYVAKGLTDSIHRYHS-DNSTKAAWDSIQSFL>
# change to 151 (original 153) and 355 (original 357)
# 2. May have multiple of the same chains (repeats) in a PDB structure and MODELLER will pick the first occurrence in the file (and I believe PDB SeqRes takes the longest one), so may need to change alignments
# 3. model's molpdf was possibly nan or > 1e8
# 4. if have a model with missing coordinates, need to rerun using a modified build_model_change_optimization function in run_modeller_uncompleted.py
# 5. can use the build_model function in run_modeller_uncompleted.py to figure out the issue
# The alignment files (.ali) are located in the the 'modeller/alignments' dir
# Once you've fixed the issues in the .ali files, rerun MODELLER as follows:
python3 run_modeller_uncompleted.py -s <script_dir>

# ----------PREDICT EDGOTYPES OF MISSENSE MUTATIONS BASED ON MUTATION LOCATION, FOLDX DDG, AND RSA----------
# NOTE: the following section should still be run on your server within your modeller script_dir (/home/username/projects/def-*/username/modeller)

# **1. Download FoldX **
# Sign up for an account on the FoldX website (https://foldxsuite.crg.eu/) to get a license 
# Once you've received a license, download the FoldX zipped tar file (e.g. _foldx5Linux64.tar_.gz)
# There should be an executable file called foldx_202x1231 (x is the year, e.g. foldx_20231231) within the zipped file

# **2. Copy the FoldX zipped tar file (e.g. _foldx5Linux64.tar_.gz) to the following directories:**
# 1. /home/username/scratch/foldx_pssm_all
# 2. /home/username/scratch/foldx_buildmodel_all

# **3. Prepare and process the output .pdb files from MODELLER for input into the FoldX PSSM command**
sbatch prepare_modeller_pdb_files_for_foldx.slurm # change the -s (--script_dir), -c (--scratch_dir), -a (--account), -f (--foldx_executable_name), and sbatch --mail_user and --account arguments

# **4. Run FoldX PSSM to determine mutation induced change in binding free energies**
cd /home/username/scratch/foldx_pssm_all # replace username 
for i in {0..n} # n is the last split_* file
do
	cp _foldx5Linux64.tar_.gz split_$i
	echo copied compressed tar to split_$i
	tar -xf split_$i/_foldx5Linux64.tar_.gz -C split_$i
	echo extracted files from compressed tar to split_$i
	cd split_$i
	sbatch foldx_pssm_$i.slurm
	echo submitted .slurm file for split_$i
	cd ../
done

# **5. Determine edgetic mutations based on output of FoldX PSSM jobs**
cd /home/username/projects/def-*/username/modeller # go back to script dir
sbatch get_edgetic_mutations.slurm # change the -s (--script_dir), -c (--scratch_dir), and sbatch --mail_user and --account arguments

# **6. Remove all files within the foldx_pssm_all dir to clear up space for foldx_buildmodel_all**
# NOTE: the foldx_pssm_all dir will have hundreds of thousands of folders/files and the scratch dir has a limit on the number of folders/files you can store
rm -r /home/username/scratch/foldx_pssm_all

# **7. Run MODELLER to build additional models needed for FoldX BuildModel**
# these models are for IR mutations that ended up being non-edgetic
python3 run_modeller_additional.py -s <script_dir> -c <scratch_dir> -p <pdb_download_dir>
DIR=_insert_your_script_directory # e.g. DIR=/home/username/projects/def-*/username/modeller
rm "$DIR"/*.ini "$DIR"/*.rsr "$DIR"/*.sch "$DIR"/*.D00000001 "$DIR"/*.V99990001
mv $DIR/*.pdb "$DIR"/results

# **8. Process additional mutations (IR mutations that are non-edgetic) to run FoldX BuildModel on**
python3 get_additional_foldx_buildmodel_mutations.slurm # change the -s (--script_dir), -c (--scratch_dir), -a (--account), -f (--foldx_executable_name), and sbatch --mail_user and --account arguments

# **9. Run FoldX BuildModel to determine mutation-induced change in protein stability/folding**
cd /home/username/scratch/foldx_buildmodel_all
for i in {0..n} # n is the last split_* file
do
	cp _foldx5Linux64.tar_.gz split_$i
	echo copied compressed tar to split_$i
	tar -xf split_$i/_foldx5Linux64.tar_.gz -C split_$i
	echo extracted files from compressed tar to split_$i
	cd split_$i
	sbatch foldx_buildmodel_$i.slurm
	echo submitted .slurm file for split_$i
	cd ../
done
for i in {0..n} # n is the last additional_* file
do
	cp _foldx5Linux64.tar_.gz additional_$i
	echo copied compressed tar to additional_$i
	tar -xf additional_$i/_foldx5Linux64.tar_.gz -C additional_$i
	echo extracted files from compressed tar to additional_$i
	cd additional_$i
	sbatch foldx_buildmodel_additional_$i.slurm
	echo submitted .slurm file for additional_$i
	cd ../
done

# **10. Install DSSP**
# if you're using Anaconda, you can install it as follows:
conda install -c salilab dssp

# **11. Calculate the relative solvent accessbility (RSA) of the mutations**
# using the procedure from Sydykova el al. (10.12688/f1000research.12874.2)
cd /home/username/projects/def-*/username/modeller # go back to <script_dir>
python3 get_rsa_slurm_file.py -s <script_dir> -c <scratch_dir> -a <acount>
for i in calc_rsa_split_*.slurm
do
sbatch $i
done
for i in calc_rsa_additional_*.slurm
do
sbatch $i
done

# **12. Categorize non-edgetic mutations as quasi-null or quasi-wildtype based on folding DDG and RSA**
sbatch get_quasi_null_wildtype.slurm # change the -s (--script_dir), -c (--scratch_dir), and sbatch --mail_user and --account arguments

# **13. Copy the following files back to local data/processed/edgotypes directory**
# (1) *_mutation_edgotypes_quasi_null_wildtype.tsv (lists mutations as edgetic, quasi-null or quasi-wildtype)
# (2) *_mutation_edgotypes.tsv (lists mutations as edgetic/non-edgetic)
# then go back to your local script_dir and run the following

# **14. Categorize mutations based on their locations (IR, buried-NIR, or exposed-NIR)**
python3 get_num_interfacial_buried_exposed.py

# **15. Categorize mutations based on their edgotypes (edgetic, quasi-null, or quasi-wildtype)**
python3 get_edgotype_quasi_null_wildtype_nums.py

# **16. Reorganize and rename columns to get final files specifying both mutation locations and edgotypes**
python get_final_mutation_location_edgotype_files.py 
