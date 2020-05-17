# markerfinder
Script for annotating phylogenetic marker proteins for concatenated protein phylogenies.  

This script will search a file of proteins against a set of 30 curated Hidden Markov Models that encode ribosomal proteins and RNA polymerase subunits. HMMs have been previously described in Sunagawa et al., Nat. Methods, 2013 (https://doi.org/10.1038/nmeth.2693).

The program assumes protein IDs are provided in the Prodigal format (i.e., contigname_1, contigname_2, etc). 
The program requires HMMER3 and the SeqIO package of Biopython. 

For questions or comments please contact Frank Aylward at faylward at vt.edu

### USAGE: python markerfinder.py <directory of protein .faa files> 

### Output
markerfinder.py provides several output files:

markerfinder_out.tsv         This is a tab-delimited output file that provides the annotation results. 

markerfinder_proteins.faa  This is the protein file with all proteins with best hits to the HMMs. Proteins are re-named to incorporate their annotation.

raw_output.txt          This is the parsed raw HMMER3 output. It can be used as a reference for debugging. 

markerfinder_cogs.txt                This is a cogs-formatted file that, together with the marerfinder_proteins.faa file, can be used as input for an ETE3 species tree workflow 
(http://etetoolkit.org/documentation/ete-build/).


