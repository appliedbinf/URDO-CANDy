# READ ME for CANDY
## Overview 

CANDY provides a command line interface to create customized amplicon databases using in silico PCR with the specified primer sequences. CANDY will create database files that serve as input for SMORE'D, a k-mer based classification tool. The required input for this program is a tabular format file containing the primers and their accompanying metadata. Optionally, additional FASTA files containing annotated sequences to include in the generated database can be used. There are two types of sequences that may be added in the additional FASTA files, presence/absence sequences and characterization sequences. These sequences differ in the way they are reported by SMORE'D. SMORE'D reports the number of times a read matches a presence/absence sequence. Only the most abundant of the sample characterization categories is reported e.g. , if there are 10,000 Mycoplasma pneumoniae P1-type 1 sequences and 1,000 Mycoplasma pneumoniae P1-type 2, the sample is characterized as Mycoplasma pneumoniae P1-type 1 and only the read count for Mycoplasma pneumoniae P1-type 1 is reported.   CANDY will output a FASTA file of representative sequences for the database, a profile file which maps the representative sequences to their annotations, and a configuration file. 

### Usage (Create database)

candy.py –i <primer table> [-o <prefix for output files>] [-d <directory for intermediate files>][-g <additional presence/absence seqs>] [-a <additional characterization seqs>]  [-t <num threads>] [-m <mapping file>]
[Bracketed arguments are optional]

#### Required arguments
   -i, --input              Input primer file; tab-separated, 6-column file  

#### Optional arguments
-o, --output_prefix      Takes a string to be added to the beginning of the output files (amplicons.fasta, profile.tsv, and config.txt)
-g, --pres_abs_seqs      Takes a FASTA file of presence/absence sequences to be included in the database
-a, --charac_seqs        Takes a FASTA file of sample characterization sequences to be included in the database
-d, --intermediate_dir   Is the name for the directory which holds intermediate files produced. Default name is date-time
-t, --threads            Takes an integer argument which tells CANDY how many threads to use
-m, --mapping            Takes a tab-separated, 2-column annotation mapping file

-l, --log                Takes a name for the log


### Alternative usage (update annotations in existing database)

candy.py --update -p <profile file> -m <mapping file>

#### Required arguments
  
--update			Tells CANDY to run in update mode
-p, --profile		The existing profile file to be updated
-m, --mapping 		The updated annotation mapping file 

