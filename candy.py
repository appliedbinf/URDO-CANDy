#!/bin/env python3

import getopt
import sys
import os
import shutil
import re
import subprocess
from multiprocessing import Pool
import logging
import json
from time import sleep


###############################################################################################
__createDB__ = True
__update_annotations__ = False
INPUT_FILE = None
__use_output_prefix__ = False
OUTPUT_PREFIX = None
__make_intermediate_dir__ = False
OUTPUT_DIR = None
__add_seqs__ = False
ADD_SEQS_FILES = None
__add_alleles__ = False
ALLELES_FILE = None
__virus_mapping__ = False
MAPPING_FILE = None
PROFILE_FILE = None
THREADS = 1
__log__=False
LOG = ''
__dependency_check__= True
__profile__ = False
# global primer_dict
primer_dict = {}
org_dict = {}



def read_primer_file(INPUT_FILE):
	"""
	Input: primer file, 6-column tab-delimited, [1. Primer name, 2. Target organism taxid(can be more than one but must be comma sep, 3. Target organism (never actually used by script), 4. Virus/Bacteria category, 5. Forward primer seqeunce, 6. reverse primer sequence]
	Output: primer_dict , org_dict
	Description: Reads the primer file and creates global dictionaries
	Checks: 1. that the input primer file has 6 columns 2. There are no repeat primer names in column 3. Column 2 must have numbers 4. Column 4 has either bacteria or virus
	"""
	logging.debug("\tPROCESS: Reading input file and generating dictionaries.")
	primers = open(INPUT_FILE, 'r').read()
	entries = [x for x in primers.split('\n') if len(x) != 0]
	#Creates a primer dictionary
	primer_dict_keys = ('primer_name', 'taxid' ,'organism','taxa', 'forward_primer_seq', 'reverse_primer_seq')
	for entry in entries:
		entry = entry.rstrip().split('\t')
		if len(entry) < 6:
			sys.exit("ERROR: Check that input primer file has the following 6 columns: \n\t1)Primer name\n\t2)Target organism taxid\n\t3)Target organism\n\t4)Virus/bacteria category\n\t5)Forward primer sequence\n\t6)Reverse complement primer sequence\n")
		entry[0] =  re.sub(pattern = ' ', string = entry[0], repl = '_')
		entry[1] = re.sub(pattern = ' ', string = entry[1], repl = '')
		if entry[0] in primer_dict:
			sys.exit(f"\nERROR: Duplicate primer name: {entry[0]}\n")
		if not bool(re.match('^[0-9,]+$',entry[1])):
			sys.exit("\nERROR: Column 2 of input primer file must be taxa ID(s) of target organism (comma separated if more than one ID.)\n")
		if not (bool(re.match("Bacteria",entry[3], flags=re.I)) or bool(re.match("Virus",entry[3], flags=re.I))):
			sys.exit("\nERROR: Column 4 must indicate the target category as either \"virus\" or \"bacteria\".\n")
		primer_dict[entry[0]] = dict(zip(primer_dict_keys, entry))
	#Create and organism dictionary
	for key in primer_dict:
		org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
		if primer_dict[key]['organism'] not in org_dict:
			org_dict[primer_dict[key]['organism']] = {}
			org_dict[primer_dict[key]['organism']]['taxa'] = primer_dict[key]['taxa']
			org_dict[primer_dict[key]['organism']]['org_name'] = org_name
			org_dict[primer_dict[key]['organism']]['taxid'] = primer_dict[key]['taxid']


def make_primer_file(primer_dict):
	"""
	Input: primer_dict
	Output: FASTA file for each primer pair [*_primer.fasta]
	Description: reads primer dict and writes primers to FASTA files as forward and reverse complement of reverse primer
	Checks: None
	"""
	#for each primer create a fasta file with the forward and reverse compliment primers
	logging.debug("\tPROCESS: Creating primer files.")
	for key in primer_dict:
		primer_name = org_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
		file_out = "{}_primers.fasta".format(primer_name)
		if os.path.exists(file_out) and os.path.isfile(file_out):
			os.remove(file_out)
		with open(file_out, 'a') as out_handle:
			logging.debug(f"\tPROCESS: Creating primer file for:{primer_name} ")
			#print(primer_dict[key]['primer_name']," Reverse: ", primer_dict[key]['reverse_primer_seq'])
			primer_dict[key]['reverse_primer_seq'] = rev_comp(primer_dict[key]['reverse_primer_seq'])
			#print(primer_dict[key]['reverse_primer_seq'])
			out_handle.write(">{primer_name}_F\n{forward_primer_seq}\n>{primer_name}_R\n{reverse_primer_seq}\n".format(**primer_dict[key]))


def create_taxid_list(org_dict):
	"""
	Input: org_dict
	Output: List of taxids [*_taxid.txt] 
	Description: Uses **taxonkit** to generate a list of all children taxids for the provided taxids in primer file 
	Checks: None
	"""
	logging.debug("\tPROCESS: Creating lists of taxa IDs")
	for key in org_dict:
		org_name = org_dict[key]['org_name']
		if os.path.isfile(f"{org_dict[key]['taxid'].replace(',','.')}_taxid.txt"):
			logging.debug(f"Taxid list already downloaded for {org_dict[key]['taxid']}, moving to next organism.")
		else:
			try:
				get_taxid_cmd = f"taxonkit list --ids {org_dict[key]['taxid']} --indent \"\" -o {org_dict[key]['taxid'].replace(',','.')}_taxid.txt"
				logging.debug(f"\tPROCESS: Creating list of taxa ids for: {org_name}")
				subprocess.call([get_taxid_cmd], shell=True)
			except subprocess.CalledProcessError as TE:
				logging.error(f"ERROR: Taxonkit unable to run with error: {TE}")
				sys.exit()
			taxid_cmd2 = f"sed -i '/^$/d' {org_dict[key]['taxid'].replace(',','.')}_taxid.txt"
			subprocess.call([taxid_cmd2], shell=True)

def blast_primers(key):
	"""
	Input: primer_dict keys [*_primer.fasta, *_taxid.txt]
	Output: file of blast results for each primer [*_blast.results]
	Description: Takes the primer fasta files and taxid lists and BLASTs primers against taxa specific sequences in nt_v5 database, output is sorted to make parsing easier
	Checks: None
	"""
	logging.debug("\tPROCESS: Starting primer BLAST")
	org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
	primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
	if os.path.exists(f"{primer_name}_blast.results") and os.path.isfile(f"{primer_name}_blast.results"):
		os.remove(f"{primer_name}_blast.results")
	try:
		blast_cmd = f"blastn -query {primer_name}_primers.fasta -db nt_v5 -taxidlist {primer_dict[key]['taxid'].replace(',','.')}_taxid.txt -outfmt '6 qseqid sseqid qlen slen length pident mismatch gaps gapopen evalue bitscore qstart qend sstart send sstrand' -word_size 7 -evalue 500000 -max_target_seqs 20000 -num_threads 10 -out {primer_name}_blast.results"
		logging.debug(f"\tPROCESS: {blast_cmd}")
		blast_pipes = subprocess.Popen(blast_cmd, shell=True,  stderr=subprocess.PIPE)
		std_err = blast_pipes.communicate()
	except subprocess.CalledProcessError as subprocess_error:
		logging.error()

	sort_cmd = f"sort -k14n {primer_name}_blast.results -o {primer_name}_blast.results"
	subprocess.run(sort_cmd, shell=True)


def extract_subsequence(genome, start, stop):
	"""
	Input: genome, start and stop of amplicons identified in BLAST parsing
	Output: returns FASTA sequences to parse_blast_results
	Descripton: Takes the gi accession number from Blast results and start and stop position and used blastdbcmd to extract seqeunce from nt_v5 database then replaces header with NCBI taxonomy using JGI http: query
	"""
	thisStart = min(start, stop)
	thisStop  = max(start, stop)

	extract_cmd = f"blastdbcmd -db nt_v5 -entry {genome} -range {thisStart}-{thisStop}"
	extract_seq = subprocess.check_output(extract_cmd, shell=True,universal_newlines=True)
	header,seq = extract_seq.split('\n',1)
	seq = seq.rstrip()
	if thisStart != start :
		seq = rev_comp(seq)

	taxonomy = re.split(">|:", header)[1]
	taxonomy_cmd = f'curl -s http://taxonomy.jgi-psf.org/sc/simple/header/{taxonomy}'
	header = ">"+ taxonomy +"|" + subprocess.check_output(taxonomy_cmd, shell=True,universal_newlines=True)
	seq = re.sub(pattern="\n", string = seq, repl="")
	return("{}\n{}\n".format(header,seq))


def parse_blast_results(key):
	"""
	Input: blast results [*_blast.results]
	Output: extracted amplicons [*_amplicons.fasta]
	Description: This fucntion parses the blast results and creates a results dictionary with the key being the matched genome which holds a list of forward and reverse coordinates, finds all forward and reverse pairs that fall within the parameters and uses the extract_seq command to pull the region from nt_v5. 
	Checks: None
	"""
	#for each blast results file create a results dictionary, each entry is a genome hit
	count = 0
	# results will be a 4-dimensional variable that will store sorted subject genomic positions as follows
	# Level 1: list - which genome the hit was found
	# Level 2: list - which primer was it - forward or reverse
	# Level 3 and 4: 2-d list - each row is the genomic match coordinates, first column is start and second column is the stop
	results = {}

	primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
	org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
	blast_file_name = primer_name + "_blast.results"
	amplicon_file = primer_name + "_amplicons.fasta"

	logging.debug(f"\tPROCESS: Creating amplicon file for: {primer_name}")

	if os.path.exists(amplicon_file) and os.path.isfile(amplicon_file):
		os.remove(amplicon_file)

	with open(amplicon_file, 'w') as out_handle:
		#print("Opening amplicon file:"+amplicon_file)
		# read the blast file
		with open(blast_file_name, 'r') as in_handle:
			for line in in_handle:
				col = line.rstrip().split('\t')
				#col[1] is a genome in the database (subject)
				col[1] = col[1].split("|")[1]

				if col[1] not in results:
					results[col[1]] = {}

				# inspect only those primers whose alignment covers >= 80% length of that primer
				if (  ( 1+abs(int(col[12]) - int(col[11])) ) / int(col[2])  ) >= .7:

					#col[0] is the primer which matched genome in db (query)
					if col[0].endswith('_F'):
						if 'F' not in results[col[1]]:
							results[col[1]]['F'] = list()

						results[col[1]]['F'].append([int(col[13]), int(col[14]), col[15]])
					elif col[0].endswith('_R'):
						if 'R' not in results[col[1]]:
							results[col[1]]['R'] = list()

						results[col[1]]['R'].append([int(col[13]), int(col[14]), col[15]])
		# read the BLAST results, processing them now
		for genome in results:
			rev_index = 0
			fwd_index = 0
			if 'F' not in results[genome] or 'R' not in results[genome]:
				continue

			while fwd_index < len(results[genome]['F']) and rev_index < len(results[genome]['R']):
				minF = min(results[genome]['F'][fwd_index][0], results[genome]['F'][fwd_index][1])

				for current_rev in range(rev_index, len(results[genome]['R'])):

					maxR = max(results[genome]['R'][current_rev][0],results[genome]['R'][current_rev][1])

						#Print statements for debugging
					# print("Range: ",range(rev_index, len(results[genome]['R'])))
					# print("fwd:", fwd_index, ", current_rev:", current_rev, ", rev_index:", rev_index, ", lenF:", len(results[genome]['F']), ", lenR:",len(results[genome]['R']))
					# print("genome ["+genome+"]. fwd primer coord:",results[genome]['F'][fwd_index],"",results[genome]['R'][current_rev],abs(maxR - minF))

					# if maxR - minF < -300 then the reverse primer binds way before forward.
					# skip further comparisons and increment the rev_index
					if maxR - minF < -300 :
						rev_index += 1
						# Go to next rev_index in range()
						continue

					# if maxR - minF > 300 then the reverse is way ahead of forward.
					# skip further comparisons and increment the fwd_index
					elif maxR - minF > 300 :
						#fwd_index += 1
						continue

					# now look if the orientations match
					# first equality is for orientation
					# second condition is to ensure the minimum primer length is 75bp
					elif results[genome]['F'][fwd_index][2] == results[genome]['R'][current_rev][2] and abs(maxR - minF) > 75:
							#print("Found overlap: ","genome ["+genome+"]. fwd primer coord:",results[genome]['F'][fwd_index],"",results[genome]['R'][current_rev],abs(maxR - minF))

						if results[genome]['F'][fwd_index][2] == "minus":
							maxR = min(results[genome]['R'][current_rev][0],results[genome]['R'][current_rev][1])
							minF = max(results[genome]['F'][fwd_index][0], results[genome]['F'][fwd_index][1])

						if  re.match("bacteria",primer_dict[key]['taxa'], flags=re.I):
							# bacteria
							out_handle.write(extract_subsequence(genome, minF, maxR))
						elif  re.match("virus",primer_dict[key]['taxa'], flags=re.I):
							# virus
							out_handle.write(extract_subsequence(genome, minF, maxR))
					else :
						# the only time you will be here is when your primers are either on different strand or too close to each other
						# print("Strand mismatch.  You shouldn't see this message repeated consecutively multiple times.")
						continue
						
				# if no hits are found for the current fwd_index, then increment it
				fwd_index += 1

	out_handle.close()

	
def rev_comp(dna):
	complement = {'A':'T','C':'G','G':'C','T':'A','W':'S','S':'W','R':'Y','Y':'R','M':'K','K':'M','N':'N','V':'V','H':'H','D':'D','B':'B',"\n":""}
	return ''.join([complement[base] for base in dna[::-1]])

		
def derep_amplicons(primer_dict):
	"""
	Input: amplicons files [*_amplicons.fasta]
	Output: dereplicated files [*_derep_amplicons.fasta]
	Description: Uses VSEARCH to dereplicate (remove all identical seqeucnes) all amplicon files.
	Checks:
	"""
	for key in primer_dict:
		primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
		amplicon_file = primer_name + "_amplicons.fasta"
		derep_file = primer_name + "_derep_amplicons.fasta"
		# comment the next line when you are actually running the program.  This is placed here for debugging purposes only
		if os.path.exists(derep_file) and os.path.isfile(derep_file):
			os.remove(derep_file)

		logging.debug(f"\tPROCESS: Dereplicating amplicons for: {primer_name}")
		dereplicate_cmd = f"vsearch --derep_fulllength {amplicon_file} --strand both --fasta_width 0 --notrunclabels --output {derep_file}".split(" ")
		subprocess.run(dereplicate_cmd, stderr=subprocess.DEVNULL)
		
		
def fix_headers(pimer_dict):
	"""
	Input: dereplicated amplicon files, (opt) mapping file [*_derep_amplicons.fasta]
	Output: combined taxonomy file [all_with_taxonomy.fasta -> taxonomy_derep.fasta]
	Description: This file parses the full NCBI taxonomy from the dereplicated amplicon files and pares it down to just species level and below then feeds everything into a single file (all_with_taxonomy.fasta. If a mapping file is provided, replaces any species level classification with the desired annotation provided in the mapping file. This fuction also dereplicates the all_with_taxonomy.fasta file to taxonomy_derep.fasta) 
	Checks: None
	"""
	logging.debug("\tPROCESS: Relabeling representative seqeunces and combing amplicons into one file")
	taxonomy_check_file = "all_with_taxonomy.fasta"
	if __virus_mapping__:
		mapping_file = MAPPING_FILE
		mapping_dict = {}
		with open(mapping_file, 'r') as in1_file:
			for line in in1_file:
				mapping_dict[line.split('\t')[0]] = line.split('\t')[1]
		with open(taxonomy_check_file, 'w') as out_file:
			for key in primer_dict:
				primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
				org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
				derep_file = primer_name + "_derep_amplicons.fasta"
				with open(derep_file, 'r') as in_file:                                
					for line in in_file:
						if line.startswith(">"):
							#print(line.rstrip())
							if "Not found" in line:
								continue
							line = line.split("s:",1)[1].rstrip()
							#print(line)
							if line in mapping_dict.keys():
								#print(">" + mapping_dict[line].rstrip())
								out_file.write(">" + mapping_dict[line])
							else:
								#print(">" + line.rstrip())
								out_file.write(">" + line + "\n")
						else:
							#print(line)
							out_file.write(line)
	else:
		with open(taxonomy_check_file, 'w') as out_file:
			for key in primer_dict:
				primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
				org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
				derep_file = primer_name + "_derep_amplicons.fasta"
				with open(derep_file, 'r') as in_file:                                
					for line in in_file:
						if line.startswith(">"):
							if "Not found" in line:
									continue
							line = line.split("s:")[1]
							out_file.write(">" + line)
						else:
							out_file.write(line)
	out_file.close()
	derep_all_cmd = f"vsearch --derep_fulllength {taxonomy_check_file} --output taxonomy_derep.fasta --fasta_width 0 --strand both --notrunclabels --minseqlength 80".split(" ")
	logging.debug("\tPROCESS: Dereplicating combined amplicon file")
	subprocess.run(derep_all_cmd, stderr=subprocess.DEVNULL)

def create_smored_files(primer_dict):
	"""
	Input: dereplicated combined amplicon file, (opt)pres/abs seqs, *(opt)charac seqs 
	Output: profile file, amplicons file, config file
	Description: This program takes the complete dereplicated collection of amplicons created with in silico PCR and combines them with any addition amplicons provided in fasta files, relabels the seqeunces with amplicon_N or their type and creates the profile file that maps the seqs to their annotation and creates a config file.  Also removes and in silico PCR generated sequencs that have ambiguous bases.
	Checks: None
	"""
	#to remove seqeunces with n reads in all_with_taxonomy.fasta
	logging.debug("\tPROCESS: Generating SMORE'D database files.")
	taxonomy_check_file = "taxonomy_derep.fasta"
	taxonomy_file = "taxonomy_ready.fasta"
	fix_taxonomy_check_cmd = f"sed ':a;N;$!ba;s/\n[^>]/PLACEHOLDER/g' {taxonomy_check_file} |grep -v 'nnn' |sed 's/PLACEHOLDER/\n/' > {taxonomy_file}"
	amplicon_dict = {}
	if __use_output_prefix__ :
		amplicon_file = OUTPUT_PREFIX + "_amplicons.fasta"
		profile_file = OUTPUT_PREFIX + "_profile.tsv"
		config_file = OUTPUT_PREFIX + "_config.txt"
	else:
		amplicon_file = "amplicons.fasta"
		profile_file = "profile.tsv"
		config_file = "config.txt"
	
	with open(amplicon_file, 'w') as out_handle:
		with open(taxonomy_check_file, 'r') as in_file:
			count = 0
			for line in in_file:
				if line.startswith(">"):
					count += 1
					new_header = "amplicon_" + str(count)
					out_handle.write(">" + new_header + "\n")
					#print(">" + new_header)
					amplicon_dict[re.sub(pattern="_", string = new_header, repl = "\t")] = re.sub(pattern = ">", string = line, repl="").rstrip()
				else:
					out_handle.write(line)
					#print(line.rstrip())
				# final_gm_count = count 
		in_file.close()
		if __add_seqs__ :
			add_seqs_file = ADD_SEQS_FILES            
			with open(add_seqs_file, 'r') as in_file2:
				for line in in_file2:
					if line.startswith(">"):
						count += 1
						new_header = "amplicon_" + str(count)
						out_handle.write(">" + new_header + "\n")
						#print(">" + new_header)
						amplicon_dict[re.sub(pattern="_", string = new_header, repl = "\t")] = re.sub(pattern = ">", string = line, repl="").rstrip()
					else:
						out_handle.write(line)
			in_file2.close()
		if __add_alleles__:
			add_alleles_file = ALLELES_FILE
			with open(add_alleles_file, 'r') as in_file3:
				for line in in_file3:
					if line.startswith(">"):
						gene = re.sub(pattern=">", string = line.split('_')[0],repl="").rstrip()
						number = line.split("_")[1].rstrip()
						annotation = line.split("_")[2].rstrip()
						new_header = gene + "_" + number
						# print(">" + new_header)
						out_handle.write(">" + new_header + "\n")
						amplicon_dict[re.sub(pattern="_", string = new_header, repl = "\t")] = annotation 
					else:
						out_handle.write(line)
			in_file3.close()
		
	with open(profile_file, 'w') as out2_handle:
		for key, value in amplicon_dict.items():
			out2_handle.write('%s\t%s\n' % (key,value))
	cwd = os.getcwd()
	with open(config_file, 'w') as out3_handle:
		out3_handle.write("[loci]\namplicon\t" +cwd +"/" + amplicon_file +"\n[profile]\nprofile\t" +cwd+ "/" +profile_file+"\n")

def update_annotations():
	"""
	Input: profile file, mapping file 
	Output: an updated profile file
	Description: this function runs only with candy.py --update and takes a previously made profile file and mapping file 
	Checks:
	"""
	profile_file = PROFILE_FILE
	mapping_file = MAPPING_FILE
	mapping_dict = {}
	temp_file = "temp"
	with open(mapping_file, 'r') as in1_file:
		for line in in1_file:
			mapping_dict[line.split('\t')[0]] = line.split('\t')[1]
	with open(temp_file, 'w') as out_file:
		with open(profile_file, 'r') as in2_file:
			for line in in2_file:
				taxonomy = line.split('\t')[2].rstrip()
				if taxonomy in mapping_dict.keys():
					line = line.split('\t')[0].rstrip()+ '\t' + line.split('\t')[1].rstrip()+ '\t' + mapping_dict[taxonomy]
					out_file.write(line)
				else:
					out_file.write(line)
	os.rename(temp_file, PROFILE_FILE)
	
	
	
def clean_up(primer_dict, org_dict):
	"""
	Input: all of the intermediate files produced
	Output:  output directory with all of the intermediate files
	Description: This function moves all of the files created in the process of making database files into a directory. If a name for the directory isn't provided it gives is a date time name
	Checks:
	"""
	if not __make_intermediate_dir__ :
		global OUTPUT_DIR
		OUTPUT_DIR = subprocess.check_output('date "+%Y%m%d_%H%M"',shell=True).decode('utf-8').rstrip()
		if os.path.isdir(OUTPUT_DIR):
			shutil.rmtree(OUTPUT_DIR)
		os.mkdir(OUTPUT_DIR)
	else:
		if os.path.isdir(OUTPUT_DIR):
			shutil.rmtree(OUTPUT_DIR)
		if os.path.isfile(OUTPUT_DIR):
			os.remove(OUTPUT_DIR)
		os.mkdir(OUTPUT_DIR)
	logging.debug(f"\tPROCESS: Moving accessory file to {OUTPUT_DIR}")
	for key in primer_dict:
		primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
		try:
			shutil.move(f"{primer_name}_primers.fasta", OUTPUT_DIR)
		except:
			pass
		try:
			shutil.move(f"{primer_name}_blast.results", OUTPUT_DIR)
		except:
			pass
		try:
			shutil.move(f"{primer_name}_amplicons.fasta", OUTPUT_DIR)
		except:
			pass
		try:
			shutil.move(f"{primer_name}_derep_amplicons.fasta", OUTPUT_DIR)
		except:
			pass
	for key in org_dict:
		try:
			shutil.move(f"{org_dict[key]['taxid'].replace(',','.')}_taxid.txt", OUTPUT_DIR)
		except:
			pass
	try:
		shutil.move("all_with_taxonomy.fasta", OUTPUT_DIR)
	except:
		pass
	try:
		shutil.move("taxonomy_derep.fasta", OUTPUT_DIR)
	except:
		pass 
	
def create_log_file(LOG):
	"""
	Input: nuttin
	Output: a log file 
	Description: This just creates a log file. If a name isn't provided with a link then the file is jsut named with a date and time

	"""
	if __log__:
		logging.basicConfig(filename=f"{LOG}.log", level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
		sys.stderr.write(f"\nWriting log file to: {LOG}\n")
	else:
		LOG=subprocess.check_output('date "+%Y%m%d_%H%M"',shell=True).decode('utf-8').rstrip() +'.log'
		logging.basicConfig(filename=LOG, level=logging.DEBUG, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
		sys.stderr.write(f"\nWriting log file to: {LOG}\n")

def check_arguments():
	if __update_annotations__:
		if not __virus_mapping__ and not __profile__:
			print("ERROR: candy.py --update requires -m <mapping file> and -p <profile file>")    
			sys.exit()
			if not os.path.isfile(MAPPING_FILE):
				#print(HELP)
				print(f"\nERROR: Mapping file ({MAPPING_FILE}) does not exist.\n")
				sys.exit()
			if not os.path.isfile(PROFILE_FILE):
				#print(HELP)
				print(f"\nERROR: Profile file ({PROFILE_FILE}) does not exit.\n")
				sys.exit()
	if __createDB__:
		if not INPUT_FILE:
			print("\nERROR: Please specify input primer file.\n")
			sys.exit()
		if not os.path.isfile(INPUT_FILE):
			print(f"\nERROR: Input file ({INPUT_FILE}) does not exist.\n")
			sys.exit()
	    
		if __add_seqs__:
			if not os.path.isfile(ADD_SEQS_FILES):
				#print(HELP)
				print(f"\nERROR: Presence/absence seqs file ({ADD_SEQS_FILES}) does not exist. \n")
				sys.exit()
		if __add_alleles__:
			if not os.path.isfile(ALLELES_FILE):
				#print(HELP)
				print(f"\nERROR: Characterization seqs file ({ADD_SEQS_FILES}) does not exist.\n")
				sys.exit()
		if __virus_mapping__:
			if not os.path.isfile(MAPPING_FILE):
				#print(HELP)
				print(f"\nERROR: Mapping file ({MAPPING_FILE}) does not exist.\n")
				sys.exit()

def check_dependencies():
	"""
	Input:
	Output:
	Description:
	"""
	if __dependency_check__:
		if __createDB__:
			print("\nChecking dependecies...\n")
			devnull = open(os.devnull)
			##Checks that VSEARCH is installed
			try:
				subprocess.run(["vsearch"], stdout=devnull, stderr=devnull)
				print("VSEARCH: Passed")
			except:
				sys.stderr.write("\n WARNING: Cannot find VSEARCH. Check that VSEARCH is downloaded and in your PATH\n\n")
				sys.exit()
			
			##Checks that taxonkit will run
			try:
				#check_taxonkit_cmd = "taxonkit list --ids 85755"
				#subprocess.Popen([check_taxonkit_cmd], stdout=devnull, stderr=devnull).communicate()
				taxonkit_proc = subprocess.Popen(["taxonkit", "list", "--ids", "85755"], stdout=devnull, stderr=subprocess.PIPE)
				taxonkit_stderr = taxonkit_proc.communicate()
				print("TaxonKit: Passed")
				if taxonkit_proc.returncode != 0:
					sys.stderr.write(f"\n {taxonkit_stderr.decode('utf8')}")
					sys.exit()
			except KeyboardInterrupt:
				sys.exit()
			except subprocess.CalledProcessError as CE: 
				sys.stderr.write(f"\nPlease check your taxonkit installation\n{CE}")
				sys.exit()
			'''
			Checks that the correct version of BLAST is installed
			'''
			blast_version = str(subprocess.check_output('blastn -version',shell=True).decode('utf-8').rstrip().split('\n')[0].split(' ')[1])
			if blast_version < '2.8.1+':
				blast_path = os.path.dirname(subprocess.check_output('which blastn',shell=True).decode('utf-8').rstrip())
				print("WARNING: You must have NCBI BLAST+ version 2.8.1 installed and in your PATH")
				print(f"WARNING: You have BLAST version: {blast_version}")
				print(f"If BLAST+ 2.8.1 is in your PATH, make sure it appears before {blast_path}")
				sys.exit()
			else:
				print("BLAST version: Passed")
			##Checks that the correct version of the BLAST database is installed
			try:
				blast_db_check_cmd = "blastdbcmd -info -db nt_v5"
				subprocess.check_output(blast_db_check_cmd, shell=True)
				print("BLAST database nt_v5: Passed")
			except subprocess.CalledProcessError as blastDBerror:
				sys.stderr.write("Make sure nt_v5 is downloaded and the ennvironmental variable 'BLASTDB' is set to the directory containing nt_v5.\n")
				exit()
############################################################################################
#Program execution starts here
HELP = "To create a new database: \n\ncandy.py -i <primer table> [-o <prefix for output files>] [-d <directory for intermediate file>] [-g <additional presence/absence seqs>] [-a <additional characterization seqs>] [-t num threads] [-m <mapping file>] [-l <log file>]  \n\nTo update an existing profile file: \n\ncandy.py --update -p <profile file> -m <mapping file> \n\nUse candy -h  for more infomation\n\n"

HELP_LONG = "To create a new database: \ncandy.py -i <primer table> [-o <prefix for output files>] [-d <directory for intermediate file>] [-g <additional presence/absence seqs>] [-a <additional characterization seqs>] [-t num threads] [-m <mapping file>] [-l <log file>]  \n\nTo update an existing profile file: \ncandy.py --update -p <profile file> -m <mapping file> \n\n-------------------------------------------------------------------------------------------------------------\n\nUsage (Create new database is the default for candy.py.)\n\tcandy.py -i <primer table> [-o <prefix for output files>] [-d <directory for intermediate file>] [-g <additional presence/absence seqs>] [-a <additional characterization seqs>] [-t num threads] [-m <mapping file>] [-l <log file>]\n\nRequired arguments: \n-i, --input\n  Input primer file; tab-delimited, 6-column file. Where, \n\tColumn 1: Primer name\n\tColumn 2: Taxonomy ID of target organism (comma separated if more than one)\n\tColumn 3: Target organism\n\tColumn 4: Virus/Bacteria category\n\tColumn 5: Forward primer sequence\n\tColumn 6: Reverse complement primer seqeunce\n\nOptional arguments:\n-o, --output_prefix\n  Takes a string to be added to the beginning of the output files (amplicons.fasta, profile.tsv, and config.txt)\n-g, --pres_abs_seqs\n  Takes a FASTA file of presence/absence sequences to be included in the database\n-a, --charac_seqs\n  Takes a FASTA file of sample characterization sequences to be included in the database\n-d, --intermediate_dir\n  Is the name for the directory which holds intermediate files produced. Default name is date-time\n-t, --threads\n  Takes an integer argument which tell CANDy how many threads to use\n-m, --mapping\n  Takes a tab-delimited, 2-column annotation mapping file\n\n\n-------------------------------------------------------------------------------------------------------------\n\nAlternative usage (update annotations in existing database)\n\tcandy.py --update -p <profile file> -m <mapping file>\n\nRequired arguments:\n--update\n  Tells candy to run in update mode\n-p, --profile\n  The existing profile file to be updated\n-m, --mapping\n  The updated annotation mapping file\n\n "
try:
	sys.argv[1]
except IndexError:
	print(HELP)
	sys.exit(0)
#######Checks for write permission in the directory, exits if user does not have right permissions
	try:
		testfile = tempfile.TemporaryFile(dir = cwd)
		testfile.close()
	except IOError:
		sys.exit('You do not have write permission in this directory')
#Input arguments
__options__, __remainders__ = getopt.getopt(sys.argv[1:], 'i:o:d:g:t:a:p:m:l:hn',[    'input=',
	'output_prefix=',
	'intermediate_directory=',
	'pres_abs_seq=',
	'threads=',
	'charac_seqs=',
	'update',
	'create',
	'profile=',
	'mapping=',
	'log=',
	'help',
	'depen_check'])

for opt, arg in __options__:
	if opt in ('-h', '-help'):
		print(HELP_LONG)
		exit()
	if opt in ('-i', '--input'):
		INPUT_FILE = arg 
	elif opt in ('-o','--output_prefix'):
		__use_output_prefix__ = True
		OUTPUT_PREFIX = arg
	elif opt in ('-d', '--intermediate_directory'):
		__make_intermediate_dir__ = True
		OUTPUT_DIR = arg
	elif opt in ('-g', '--pres_abs_seqs'):
		__add_seqs__ = True
		ADD_SEQS_FILES = arg
	elif opt in ('-t', '--threads'):
		try:
			THREADS = int(arg)
		except ValueError:
			print("Error: Enter an interger value for threads.")
			sys.exit(0)
	elif opt in ('-a', '--charac_seqs'):
		__add_alleles__= True
		ALLELES_FILE = arg
	elif opt in ('--update'):
		__update_annotations__ = True
		__createDB__ = False
	elif opt in ('-p', '--profile'):
		__profile__ = True
		PROFILE_FILE = arg
	elif opt in ('-m', '--mapping'):
		__virus_mapping__ = True
		MAPPING_FILE = arg
	elif opt in ('-l','--log'):
		LOG = arg
		__log__= True
	elif opt in ('-n','--depen_check'):
		__dependency_check__ = False        

import tempfile
import errno    
def main():
	check_arguments()
	check_dependencies()
	create_log_file(LOG)
	if __createDB__ :     
		read_primer_file(INPUT_FILE)
		#make_primer_file(primer_dict)
		#create_taxid_list(org_dict)
		with Pool(int(THREADS)) as pool:
			#pool.map(blast_primers, primer_dict.keys())
			pool.map(parse_blast_results, primer_dict.keys())
		derep_amplicons(primer_dict)
		fix_headers(primer_dict)
		create_smored_files(primer_dict)
		clean_up(primer_dict, org_dict)
	elif __update_annotations__:
		update_annotations()
if __name__ == '__main__':
	main()

###Checks needed
#Check user's blast version and blast database version blastdbcmd -version |head -1 |cut -d':' -f2 |sed 's/\+//' | sed 's/^\s//'
# get path to blastdatabase
#make newest blast cmd