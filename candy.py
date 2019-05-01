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

# global primer_dict
primer_dict = {}
org_dict = {}



def read_primer_file(INPUT_FILE):
    logging.debug("Analysis: Reading input file and generating dictionaries.")
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
    logging.debug("Analysis: Creating primer files.")
    for key in primer_dict:
        primer_name = org_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
        file_out = "{}_primers.fasta".format(primer_name)
        if os.path.exists(file_out) and os.path.isfile(file_out):
            os.remove(file_out)
        with open(file_out, 'a') as out_handle:
            logging.debug(f"Creating primer file for:{primer_name} ")
            out_handle.write(">{primer_name}_F\n{forward_primer_seq}\n>{primer_name}_RC\n{reverse_comp_primer_seq}\n".format(**primer_dict[key]))


def create_taxid_list(org_dict):
    logging.debug("Analysis: Creating lists of taxa IDs")
    for key in org_dict:
        org_name = org_dict[key]['org_name']
        if os.path.isfile(f"{org_dict[key]['taxid'].replace(',','.')}_taxid.txt"):
            logging.debug(f"Taxid list already downloaded for {org_dict[key]['taxid']}, moving to next organism.")
        else:
            try:
                get_taxid_cmd = f"taxonkit list --ids {org_dict[key]['taxid']} --indent \"\" -o {org_dict[key]['taxid'].replace(',','.')}_taxid.txt"
                logging.debug(f"Analysis: Creating list of taxa ids for: {org_name}")
                subprocess.call([get_taxid_cmd], shell=True)
            except subprocess.CalledProcessError as TE:
                logging.error(f"ERROR: Taxonkit unable to run with error: {TE}")

            #     taxid_pipes = subprocess.Popen(get_taxid_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #     std_out, std_err = taxid_pipes.communicate()
            #     logging.debug(f"{std_out.decode('utf-8')}{std_err.decode('utf-8')}")
            # except subprocess.CalledProcessError as subprocess_error:
            #     logging.error(f"Could not download taxid list for: {org_dict[key]['taxid']}")
            #     logging.error(f"ERROR: {subprocess_error}")
            #     sys.exit(f"Could not download taxid for: {org_dict[key]}")

            taxid_cmd2 = f"sed -i '/^$/d' {org_dict[key]['taxid'].replace(',','.')}_taxid.txt"
            subprocess.call([taxid_cmd2], shell=True)
def blast_primers(key):
    #blast primer files against genome database
    logging.debug("Analysis: Starting primer BLAST")
    org_name = re.sub(pattern = ' ', string = primer_dict[key]['organism'].rstrip(), repl='_')
    primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
    if os.path.exists(f"{primer_name}_blast.results") and os.path.isfile(f"{primer_name}_blast.results"):
        os.remove(f"{primer_name}_blast.results")
    try:
        blast_cmd = f"blastn -query {primer_name}_primers.fasta -db nt_v5 -taxidlist {primer_dict[key]['taxid'].replace(',','.')}_taxid.txt -outfmt '6 qseqid sseqid qlen slen length pident mismatch gaps gapopen evalue bitscore qstart qend sstart send sstrand' -word_size 7 -evalue 500000 -max_target_seqs 20000 -num_threads 10 -out {primer_name}_blast.results"
        logging.debug(f"Analysis: \n{blast_cmd}")
        blast_pipes = subprocess.Popen(blast_cmd, shell=True,  stderr=subprocess.PIPE)
        std_err = blast_pipes.communicate()
    except subprocess.CalledProcessError as subprocess_error:
        logging.error()

    sort_cmd = f"sort -k14n {primer_name}_blast.results -o {primer_name}_blast.results"
    subprocess.call(sort_cmd, shell=True)


def extract_subsequence(genome, start, stop):
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

    logging.debug(f"Creating amplicon file for: {primer_name}")

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
                        
                # if no hits are found for the current fwd_index, then increment it
                fwd_index += 1

    out_handle.close()

    
def rev_comp(dna):
    complement = {'A':'T','C':'G','G':'C','T':'A','W':'S','S':'W','R':'Y','Y':'R','M':'K','K':'M','N':'N','V':'V','H':'H','D':'D','B':'B',"\n":""}
    return ''.join([complement[base] for base in dna[::-1]])

        
def derep_amplicons(primer_dict):
    #this fucntion runs vsearch cluster to remove identical amplicons
    for key in primer_dict:
        primer_name = re.sub(pattern = ' ', string = primer_dict[key]['primer_name'].rstrip(), repl='_')
        amplicon_file = primer_name + "_amplicons.fasta"
        derep_file = primer_name + "_derep_amplicons.fasta"
        # comment the next line when you are actually running the program.  This is placed here for debugging purposes only
        if os.path.exists(derep_file) and os.path.isfile(derep_file):
            os.remove(derep_file)

        logging.debug(f"Dereplicating amplicons for: {primer_name}")
        dereplicate_cmd = f"vsearch --derep_fulllength {amplicon_file} --strand both --fasta_width 0 --notrunclabels --output {derep_file}".split(" ")
        subprocess.run(dereplicate_cmd, stderr=subprocess.DEVNULL)
        
        #os.system(dereplicate_cmd)

        
def fix_headers(pimer_dict):
    logging.debug("Analysis: Relabeling representative seqeunces and combing amplicons into one file")
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
                            line = line.split("s:")[1].rstrip()
                            if line in mapping_dict.keys():
                                out_file.write(">" + mapping_dict[line])
                            else:
                                out_file.write(">" + line)
                        else:
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
                            line = line.split("s:")[1]
                            out_file.write(">" + line)
                        else:
                            out_file.write(line)
    out_file.close()
    derep_all_cmd = f"vsearch --derep_fulllength {taxonomy_check_file} --output taxonomy_derep.fasta --fasta_width 0 --strand both --notrunclabels --minseqlength 80".split(" ")
    logging.debug("Analysis: Dereplicating combined amplicon file")
    subprocess.run(derep_all_cmd, stderr=subprocess.DEVNULL)
    #os.system(derep_all_cmd)

final_gm_count = 0
def create_smored_files(primer_dict):
    #to remove seqeunces with n reads in all_with_taxonomy.fasta
    logging.debug("Analysis: Generating SMORE'D database files.")
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
                    out_file.write(line)
                else:
                    out_file.write(line)
    os.rename(temp_file, PROFILE_FILE)
    
    
    
def clean_up(primer_dict, org_dict):
# add a step to clean up directory to remove BLAST output, log files and derep files
    if __make_intermediate_dir__ :
        if os.path.isdir(OUTPUT_DIR):
            shutil.rmtree(OUTPUT_DIR)
        os.mkdir(OUTPUT_DIR)
    else:
        OUTPUT_DIR = subprocess.check_output('date "+%Y%m%d_%H%M"',shell=True).decode('utf-8').rstrip()
        if os.path.isdir(OUTPUT_DIR):
            shutil.rmtree(OUTPUT_DIR)
        os.mkdir(OUTPUT_DIR)
    logging.debug(f"Moving accessory file to {OUTPUT_DIR}")
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
    
def create_log_file():
    if __log__:
        logging.basicConfig(filename=LOG, level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        sys.stderr.write(f"\nWriting log file to: {LOG}\n")
    else:
        LOG=subprocess.check_output('date "+%Y%m%d_%H%M"',shell=True).decode('utf-8').rstrip() +'.log'
        logging.basicConfig(filename=LOG, level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
        sys.stderr.write(f"\nWriting log file to: {LOG}\n")

#def check_arguments():

def check_dependencies():
    if __createDB__:
        devnull = open(os.devnull)
        ##Checks that VSEARCH is installed
        try:
            subprocess.run(["vsearch"], stdout=devnull, stderr=devnull)
        except:
            sys.stderr.write("\n WARNING: Cannot find VSEARCH. Check that VSEARCH is downloaded and in your PATH\n\n")
            sys.exit()
        
        ##Checks that taxonkit will run
        try:
            #check_taxonkit_cmd = "taxonkit list --ids 85755"
            #subprocess.Popen([check_taxonkit_cmd], stdout=devnull, stderr=devnull).communicate()
            taxonkit_proc = subprocess.Popen(["taxonkit", "list", "--ids", "85755"], stdout=devnull, stderr=subprocess.PIPE)
            taxonkit_stderr = taxonkit_proc.communicate()
            if taxonkit_proc.returncode != 0:
                sys.stderr.write(f"\n {taxonkit_stderr.decode('utf8')}")
                sys.exit()
        except KeyboardInterrupt:
            sys.exit()
        except subprocess.CalledProcessError as CE: 
            sys.stderr.write(f"\nPlease check your taxonkit installation\n{CE}")
            sys.exit()
        


        
        ##Checks that the correct version of BLAST is installed
        blast_version = str(subprocess.check_output('blastn -version',shell=True).decode('utf-8').rstrip().split('\n')[0].split(' ')[1])
        if blast_version < '2.8.1+':
            blast_path = os.path.dirname(subprocess.check_output('which blastn',shell=True).decode('utf-8').rstrip())
            print("WARNING: You must have NCBI BLAST+ version 2.8.1 installed and in your PATH")
            print(f"WARNING: You have BLAST version: {blast_version}")
            print(f"If BLAST+ 2.8.1 is in your PATH, make sure it appears before {blast_path}")
            exit()
        
        ##Checks that the correct version of the BLAST database is installed
        try:
            blast_db_check_cmd = "blastdbcmd -info -db nt_v5"
            subprocess.check_output(blast_db_check_cmd, shell=True)
        except subprocess.CalledProcessError as blastDBerror:
            sys.stderr.write("Make sure nt_v5 is downloaded and the ennvironmental variable 'BLASTDB' is set to the directory containing nt_v5.\n")
            exit()
############################################################################################
#Program execution starts here
HELP = "To create a new database: \ncandy.py -i <primer table> [-o <prefix for output files>] [-d <directory for intermediate file>] [-g <additional presence/absence seqs>] [-a <additional characterization seqs>] [-t num threads] [-m <mapping file>]  \n\nTo update an existing profile file: \ncandy.py --update -p <profile file> -m <mapping file> \n\n"
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
__options__, __remainders__ = getopt.getopt(sys.argv[1:], 'i:o:d:g:t:a:p:m:l:h',[    'input=',
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
    'help'])

for opt, arg in __options__:
    if opt in ('-h', '-help'):
        print(HELP)
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
        PROFILE_FILE = arg
    elif opt in ('-m', '--mapping'):
        __virus_mapping__ = True
        MAPPING_FILE = arg
    elif opt in ('-l','--log'):
        __log__= True
        LOG = arg
        
        

import tempfile
import errno    
def main():
    check_dependencies()
    create_log_file()
    if __createDB__ :     
        read_primer_file(INPUT_FILE)
        # make_primer_file(primer_dict)
        # create_taxid_list(org_dict)
        # with Pool(int(THREADS)) as pool:
        #     pool.map(blast_primers, primer_dict.keys())
        #     pool.map(parse_blast_results, primer_dict.keys())
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