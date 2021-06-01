# markerfinder
Script for annotating phylogenetic marker proteins for concatenated protein phylogenies.  

This script will search a file of proteins against a set of 40 curated Hidden Markov Models for phylogenetic marker genes. These markers are useful for phylogenetic analysis of both bacteria and archaea, and they have been previously described in Sunagawa et al., Nat. Methods, 2013 (https://doi.org/10.1038/nmeth.2693).

The program assumes protein IDs are provided in the Prodigal format (i.e., contigname_1, contigname_2, etc). 
The program requires HMMER3 (in your PATH) and the SeqIO package of Biopython. 

The program will search for fragmented RNA Polymerase subunits and join them in silico for concatenated alignments. These genes are sometimes fragmented in cyanobacteria, thermophilic archaea, and other microbial lineages, which can cause problems for proper annotation and incorporation into concatenated alignments. Proteins are only joined if they are encoded on the same contig/replicon and have non-overlapping HMMER alignments to the COG0085 or COG0086 HMMs. 

There are options that allow for slightly different marker gene sets to be used. 

all: all 40 marker genes

ribo: only 27 ribosomal proteins

rnap: only 3 RNAP subunits

ribo_rnap: 27 ribosomal proteins and 3 RNAP subunits


### MINIMAL USAGE: python markerfinder.py -i <directory of protein .faa files> -n project_name

### Options

**-p, --proximity**
The proximity of genes to merge is the number of genes up- and down-stream of the initial hit that the program will search for additional hits to merge. This is only done for the COG0085 and COG0086 markers (the multimeric RNAP subunits). Hits outside this range will be considered independent hits. Default is 100.

**-t, --cpus**
How many CPUs/threads to use for the hmmsearch and clustal steps

**-m, --markerset**
HMM database to use. Options are "all", "ribo", "rnap" or "ribo_rnap". See README for details

**-r, --redo**
If you have already run script and you want to re-run it with different parameters, you can use the -r flag to avoid re-running HMMER (this saves a bit of time if you're running multiple times)

**-c, --concat**
If this option is specified the script will also output a concatenated alignment of the chose marker genes (FASTA format, one entry per taxa, each marker gene aligned separately with Clustal Omega)

**-a, --allhits**
If this option is specified then all hits marker genes (that are above the predefined bit score thresholds) will be output. This option is not compatible with the -c option. This can be useful if you want to see if certain marker genes are present in multiple copies. 



### Output files
markerfinder.py provides several output files, all with the prefix designated with the -n option:

*full_output.txt         This is the main tab-delimited output file that provides the annotation results. 

*.faa  This is the protein file with all merged and unmerged proteins with best hits to the HMMs. Proteins are re-named to accommodate potential merged proteins. 

*raw_output.txt          This is the parsed raw HMMER3 output (no marker gene joining). It can be used as a reference for debugging, but you can usually ignore it. 

*.cogs.txt                This is a cogs-formatted file, in the same general format used by the ETE3 toolkit, and it's used as a reference for producing the concatenated alignment. You can also use it with the *.faa output file if you want to make a tree with the ETE3 toolkit instead (http://etetoolkit.org/documentation/ete-build/).. 

*.table.tsv              This is an occurrence table for the number of marker genes that were identified for each file in the input folder. Note that if several query proteins with hits to the same marker gene were joined they only counted once. Also, this table will only show multiple hits if the -a option is used, otherwise only best hits are recorded. 

*.concat.aln           If the -c option is chosen then a concatenated alignment is produced. This is FASTA formatted, each marker gene is aligned separately, and strings of X are used to fill spaces left by missing marker genes. This alignment is not trimmed in any way, so you may wish to process it further before phylogenetic analysis. 

log_file.txt          This is just a log file of some of the outputs produced by hmmsearch and clustal that is used for debugging purposes. 

In addition, for each .faa file in the input folder a .domout and .domout.parsed file is created. These are used if you re-run this tool with the -r flag. 


### Examples

To get the best hits to a set of ribosomal markers only using 4 threads. 
>python markerfinder.py -i test_input -n test_run -t 4 -m ribo

To get all hits, including "secondary hits", or second-best hits:
>python markerfinder.py -i test_input -n test_run -t 4 -a

To get best hits and also generate a concatenated alignment: 
>python markerfinder.py -i test_input -n test_run -t 4 -c




