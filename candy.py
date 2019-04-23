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

#primer_file = "XEEPsmored-assaylist3.txt"

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

PROFILE_FILE = None
THREADS = 1




# global primer_dict
primer_dict = {}
org_dict = {}



def read_primer_file(INPUT_FILE):
    print("#######################################################\nReading input file and generating dictionaries.\n#######################################################\n")
    primers = open(INPUT_FILE, 'r').read()
    entries = [x for x in primers.split('\n') if len(x) != 0]
    #Creates a primer dictionary
    for entry in entries:
        entry = entry.rstrip().split('\t')
        entry[0] =  re.sub(pattern = ' ', string = entry[0], repl = '_')
        primer_dict_keys = ('primer_name', 'taxid' ,'organism','taxa', 'forward_primer_seq', 'reverse_comp_primer_seq')
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
    #for each primer create a fasta file with the forward and reverse compliment primers
    print("#######################################################\nCreating primer files.\n#######################################################\n")
    for key in primer_dict:
        primer_name = org_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
        file_out = "{}_primers.fasta".format(primer_name)
        if os.path.exists(file_out) and os.path.isfile(file_out):
            os.remove(file_out)
        with open(file_out, 'a') as out_handle:
            print("Creating primer file for: ", primer_name)
            out_handle.write(">{primer_name}_F\n{forward_primer_seq}\n>{primer_name}_RC\n{reverse_comp_primer_seq}\n".format(**primer_dict[key]))


def create_taxid_list(org_dict):
    print("#######################################################\nCreating lists of taxa IDs\n#######################################################\n")
    for key in org_dict:
        org_name = org_dict[key]['org_name']
        if os.path.isfile(f"{org_dict[key]['taxid'].replace(',','.')}_taxid.txt"):
            print(f"Taxid list already downloaded for {org_dict[key]['taxid']}, moving to next organism.")
        else:
            get_taxid_cmd = f"taxonkit list --ids {org_dict[key]['taxid']} --indent \"\" | sed '/^$/d' > {org_dict[key]['taxid'].replace(',','.')}_taxid.txt"
            print(get_taxid_cmd)
            os.system(get_taxid_cmd)

            
def blast_primers(key):
    #blast primer files against genome database
    #print("#######################################################\nStarting primer BLAST\n#######################################################\n")
    #for key in primer_dict:
    org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
    primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
    if os.path.exists(f"{primer_name}_blast.results") and os.path.isfile(f"{primer_name}_blast.results"):
        os.remove(f"{primer_name}_blast.results")

    blast_cmd = f"/storage/blastdb/v5/ncbi-blast-2.8.1+/bin/blastn -query {primer_name}_primers.fasta -db /storage/blastdb/v5/nt_v5 -taxidlist {primer_dict[key]['taxid'].replace(',','.')}_taxid.txt -outfmt '6 qseqid sseqid qlen slen length pident mismatch gaps gapopen evalue bitscore qstart qend sstart send sstrand' -word_size 7 -evalue 500000 -max_target_seqs 20000 -num_threads 10| sort -k14n > {primer_name}_blast.results"
    print(blast_cmd)
    os.system(blast_cmd)


def extract_subsequence(genome, start, stop):
    thisStart = min(start, stop)
    thisStop  = max(start, stop)
    extract_cmd = f"/storage/blastdb/v5/ncbi-blast-2.8.1+/bin/blastdbcmd -db /storage/blastdb/v5/nt_v5 -entry {genome} -range {thisStart}-{thisStop}"
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
    #for each blast results file create a results dictionary, each entry is a genome hit
    #print("#######################################################\nGenerating amplicon files from BLAST results\n#######################################################\n")
    #for key in primer_dict:
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

    print("Creating amplicon file for: " + primer_name)

    if os.path.exists(amplicon_file) and os.path.isfile(amplicon_file):
        os.remove(amplicon_file)

    with open(amplicon_file, 'w') as out_handle:
        print("Opening amplicon file:"+amplicon_file)
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
                    elif col[0].endswith('_RC'):
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
                        #pass
                    # move the appropriate pointers if you are close to the end so you are not stuck
                    # if rev_index == len(results[genome]['R'])-1:
                    #         fwd_index += 1
                    # if fwd_index == len(results[genome]['F'])-1:
                    #         rev_index += 1

                # if no hits are found for the current fwd_index, then increment it
                fwd_index += 1

    out_handle.close()

    
def rev_comp(dna):
    complement = {'A':'T','C':'G','G':'C','T':'A','W':'S','S':'W','R':'Y','Y':'R','M':'K','K':'M','N':'N','V':'V','H':'H','D':'D','B':'B',"\n":""}
    return ''.join([complement[base] for base in dna[::-1]])

        
def derep_amplicons(primer_dict):
    #this fucntion runs vsearch cluster to remove identical amplicons
    print("#######################################################\nDereplicating amplicon files\n#######################################################\n")
    for key in primer_dict:
        primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
        amplicon_file = primer_name + "_amplicons.fasta"
        derep_file = primer_name + "_derep_amplicons.fasta"
        # comment the next line when you are actually running the program.  This is placed here for debugging purposes only
        if os.path.exists(derep_file) and os.path.isfile(derep_file):
            os.remove(derep_file)

        print("Dereplicating amplicons for: " + primer_name)
        dereplicate_cmd = f'vsearch --derep_fulllength {amplicon_file} --strand both --fasta_width 0 --notrunclabels --output {derep_file}'
        os.system(dereplicate_cmd)

        
def fix_headers(pimer_dict):
    taxonomy_check_file = "all_with_taxonomy.fasta"
    mapping_file = "virus_mapping.txt"
    mapping_dict = {}
    with open(mapping_file, 'r') as in1_file:
        for line in in1_file:
            mapping_dict[line.split('\t')[0]] = line.split('\t')[1]
    
    print("#######################################################\nRelabeling representative seqeunces and creating mega FASTA.\n#######################################################\n")
    with open(taxonomy_check_file, 'w') as out_file:
        for key in primer_dict:
            primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
            org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
            derep_file = primer_name + "_derep_amplicons.fasta"
            #taxonomy_check_file = primer_name + "_tax.fasta"
            with open(derep_file, 'r') as in_file:
                if  re.match("bacteria",primer_dict[key]['taxa'], flags=re.I):
                    #print("------------------------------------------")
                    #print("Checking taxonomy for: " + primer_name)
                    for line in in_file:
                        if line.startswith(">"):
                            if "Not found" in line:
                               out_file.write(">" + primer_name + "\n")
                            else:
                                line = line.split("s:")[1]
                                out_file.write(">" + line)
                                
                                # if ";ss:" in line:
                                    # line = line.split("ss:")[0]
                                    # out_file.write(">" + line)
                                # if ";" in line:
                                    # line = line.split(";")[0]
                                    # out_file.write(">" + line)
                                # else:
                                    # out_file.write(">" + line)
                            # print(line.rstrip())
                        else:
                            out_file.write(line)
                            
                elif re.match("virus",primer_dict[key]['taxa'], flags=re.I):
                    for line in in_file:
                        if line.startswith(">"):
                            line = line.split("s:")[1].rstrip()
                            if line in mapping_dict.keys():
                                out_file.write(">" + mapping_dict[line])
                            else:
                                out_file.write(">" + line)
                            # out_file.write(">" + org_name + "\n")
                        else:
                            out_file.write(line)
    out_file.close()
    derep_all_cmd = f"vsearch --derep_fulllength {taxonomy_check_file} --output taxonomy_derep.fasta --fasta_width 0 --strand both --notrunclabels --minseqlength 80"
    os.system(derep_all_cmd)

final_gm_count = 0
def create_smored_files(primer_dict):
    #to remove seqeunces with n reads in all_with_taxonomy.fasta
    print("#######################################################\nGenerating SMORE'D database files.\n#######################################################\n")
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
    profile_file = PROFILE_FILE
    mapping_file = "virus_mapping.txt"
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
                    out_file.write(line.rstrip())
                else:
                    out_file.write(line.rstrip())
    os.rename(temp_file, PROFILE_FILE)
    
    
    
def clean_up(primer_dict, org_dict):
# add a step to clean up directory to remove BLAST output, log files and derep files
    print("#######################################################\nRemoving accessory files\n#######################################################\n")
    if __make_intermediate_dir__ :
        if os.path.isdir(OUTPUT_DIR):
            os.remove(OUTPUT_DIR)
        os.mkdir(OUTPUT_DIR)
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
                shutil.move(f"{org_dict[key]['taxid']}_taxid.txt", OUTPUT_DIR)
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
    else:
        for key in primer_dict:
            primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
            try:
                os.remove(f"{primer_name}_primers.fasta")
            except:
                pass
            try:
                os.remove(f"{primer_name}_blast.results")
            except:
                pass
            try:
                os.remove(f"{primer_name}_amplicons.fasta")
            except:
                pass
            try:
                os.remove(f"{primer_name}_derep_amplicons.fasta")
            except:
                pass
        for key in org_dict:
            try:
                os.remove(f"{org_dict[key]['taxid']}_taxid.txt")
            except:
                pass
        try:
            os.remove("all_with_taxonomy.fasta")
        except:
            pass
        try:
            os.remove("taxonomy_derep.fasta")
        except:
            pass
############################################################################################
#Program execution starts here

try:
    sys.argv[1]
except IndexError:
    print("To create a new database: \ncandy.py -i <input file> [-o <database files prefix>] [-d <directory for intermediate file>] [-g <sequences to add to database>] [-t num_threads] [-a <additional seqs to add>]\n\nTo update an existing profile file: \ncandy.py --update -p <profile file>\n\n")
    sys.exit(0)

#Input arguments
__options__, __remainders__ = getopt.getopt(sys.argv[1:], 'i:o:d:g:t:a:p:',[    'input=',
    'output_prefix=',
    'intermediate_directory=',
    'pres_abs_seq=',
    'threads=',
    'charac_seqs=',
    'update',
    'create',
    'profile='])

for opt, arg in __options__:
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
        THREADS = arg
    elif opt in ('-a', '--charac_seqs'):
        __add_alleles__= True
        ALLELES_FILE = arg
    elif opt in ('--update'):
        __update_annotations__ = True
        __createDB__ = False
    elif opt in ('-p', '--profile'):
        PROFILE_FILE = arg
        
    
def main():
    if __createDB__ :
        read_primer_file(INPUT_FILE)
        make_primer_file(primer_dict)
        create_taxid_list(org_dict)
        with Pool(int(THREADS)) as pool:
            pool.map(blast_primers, primer_dict.keys())
            pool.map(parse_blast_results, primer_dict.keys())
        derep_amplicons(primer_dict)
        fix_headers(primer_dict)
        create_smored_files(primer_dict)
        clean_up(primer_dict, org_dict)
    elif __update_annotations__:
        update_annotations()
if __name__ == '__main__':
    main()

