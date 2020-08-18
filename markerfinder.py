import os, sys, subprocess, argparse, re, shlex, glob, operator
import numpy as np
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import inspect

def path_to_know_finction():
	return 

hmm_db = os.path.dirname(os.path.abspath(inspect.getfile(path_to_know_finction)))+"/hmm/embl.hmm" #"hmm/RNAP.full.hmm"
#speci_db = "hmm/all.hmm"

score_dict = {"COG0012":float(210), "COG0016":float(240), "COG0018":float(340), "COG0048":float(100), "COG0049":float(120), "COG0052":float(140), "COG0080":float(90), "COG0081":float(130), "COG0085":float(200), "COG0086":float(200), "COG0087":float(120), "COG0088":float(110), "COG0090":float(180), "COG0091":float(80), "COG0092":float(120), "COG0093":float(80), "COG0094":float(110), "COG0096":float(80), "COG0097":float(100), "COG0098":float(140), "COG0099":float(120), "COG0100":float(80), "COG0102":float(100), "COG0103":float(80), "COG0124":float(320), "COG0172":float(170), "COG0184":float(60), "COG0185":float(70), "COG0186":float(80), "COG0197":float(70), "COG0200":float(60), "COG0201":float(210), "COG0202":float(80), "COG0215":float(400), "COG0256":float(70), "COG0495":float(450), "COG0522":float(80), "COG0525":float(740), "COG0533":float(300), "COG0541":float(450), "COG0552":float(220), "COG0086":float(300)}

combined_output = open("output.txt", "w")
final_proteins = []

#################################################################
############# define hmm launcher function ######################
#################################################################
def hmm_launcher(folder):
	print "Running HMMER3..."
	for files in os.listdir(folder):
		if files.endswith(".faa"):
			#print files
			input_file = os.path.join(folder, files)	
			dom_output = re.sub(".faa", ".domout", files)
			speci_dom_output = os.path.join(folder, dom_output)

			# run against the RNAP models
			cmd = "hmmsearch --cpu 16 --domtblout "+ speci_dom_output +" "+ hmm_db + " " + input_file
			#print cmd
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("hmm.out", 'w'), stderr=open("error_file.txt", 'a'))

# end

################################################################
###### Loop through and parse the checkm HMM output ############
################################################################

def hmm_parser(folder, suffix, output):
	
	record_list = []
	score_list = {}
	prot_list = []
	combined_output.write("protein\tacc\thit\tstart\tend\taln_length\tscore\tcategory\n")
	hits = []
	bit_dict = {}

	for filenames in os.listdir(folder):
		if filenames.endswith(suffix):

			acc = re.sub(suffix, "", filenames)

			f = open(folder+"/"+filenames, 'r')
			o = open(folder+"/"+filenames+".parsed", 'w')

			faa_file = re.sub(".domout", ".faa", filenames)
			protein_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(folder, faa_file), "fasta"))

			hit_dict = {}
			start_dict = {}
			end_dict = {}
			bit_dict = defaultdict(int)
			hit_type = {}
			marker_dict = {}

			for line in f.readlines():
				if line.startswith("#"):
					pass
				else:
					newline = re.sub( '\s+', '\t', line)
					list1 = newline.split('\t')
					ids = list1[0]

					hit = list1[3]

					start = list1[15]
					end   = list1[16]
					#print start, end

					score = float(list1[7])
					if score > score_dict[hit]:

						if score > bit_dict[ids]:
							hit_dict[ids] = hit
							start_dict[ids] = start
							end_dict[ids] = end
							bit_dict[ids] = score

			bit_sorted = sorted(bit_dict.items(), key=operator.itemgetter(1), reverse=True)
			output_list = []
			for item in bit_sorted:
				entry = item[0]
				output_list.append(entry +"\t"+ str(hit_dict[entry]) +"\t"+ str(start_dict[entry]) +"\t"+ str(end_dict[entry]) +"\t"+ str(bit_dict[entry]) )

			hit_profile = defaultdict(int)
			done = []
			for line in output_list:
				line1 = line.rstrip()
				tabs = line1.split("\t")
				ids = tabs[0]
				
				record = protein_dict[ids]
				record_list.append(record)

				hits.append(ids)
				cog = tabs[1]
				start = tabs[2]
				end = tabs[3]
				aln_length = str(abs(float(end) - float(start)))

				score = tabs[4]
				nr = acc +"_"+ cog

				if nr in done:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tNH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ score +"\tNH\n")
				else:
					combined_output.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ aln_length +"\t"+ score +"\tBH\n")
					o.write(ids +"\t"+ acc +"\t"+ cog +"\t"+ start +"\t"+ end +"\t"+ score +"\tBH\n")
					done.append(nr)
			o.close()

#parse speci outputs

################################################################
########## Define function for parsing HMMER3 output ###########
################################################################
def parse_domout(path_to_parsed_hmmfile, acc, protein_dict, cog_name):
	protein2dups = defaultdict(lambda:"single_besthit")
	parsed = open(path_to_parsed_hmmfile, "r")
	done = {}
	protein2coords = defaultdict(list)
	protein2align_length = {}

	main_hit = "NAN"
	rnap_hits = []
	protein2cog = defaultdict(lambda:"NA")
	protein2acc = {}
	protein2score = {}
	protein2category = {}
	protein2length = {}

	for n in parsed.readlines():
		line = n.rstrip()
		tabs = line.split("\t")
		protein = tabs[0]
		annot = tabs[2]
		if annot == cog_name:
			rnap_hits.append(protein)
			id_hit = protein +"|"+ annot
			hmm_score = float(tabs[5])
			category = tabs[6]

			if hmm_score > 1:

				start = int(tabs[3])
				end =   int(tabs[4])

				record = protein_dict[protein]
				prot_length = len(record.seq)

				nr = acc +"_"+ annot

				protein2cog[protein]    = annot
				protein2acc[protein]    = acc
				protein2score[protein]  = hmm_score
				protein2length[protein] = prot_length

				protein2coords[id_hit].append(start)
				protein2coords[id_hit].append(end)

				align_length = abs(end - start)
				protein2align_length[id_hit] = align_length

				if category == "BH" and annot == cog_name:
					main_hit = protein
					protein2dups[id_hit]

				protein2category[protein] = category

	parsed.close()
	return main_hit, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length, protein2dups


# start program
def run_program(input, project, database, cpus):

	# pick dtabase
	if database == "rnap":
		cog_set = ["COG0085", "COG0086", "COG0202"] # 3 RNAP subunits
		print("Using the RNAP marker set")
	elif database == "ribo":
		cog_set = ["COG0012", "COG0048", "COG0049", "COG0052", "COG0080", "COG0081", "COG0087", "COG0088", "COG0090", "COG0091", "COG0092", "COG0093", "COG0094", "COG0096", "COG0097", "COG0098", "COG0099", "COG0100", "COG0102", "COG0103", "COG0184", "COG0185", "COG0186", "COG0197", "COG0200", "COG0256", "COG0522"] # 27 ribosomal proteins
		print("Using the ribosomal marker set")
	elif database == "ribo_rnap":
		cog_set = ["COG0012", "COG0048", "COG0049", "COG0052", "COG0080", "COG0081", "COG0087", "COG0088", "COG0090", "COG0091", "COG0092", "COG0093", "COG0094", "COG0096", "COG0097", "COG0098", "COG0099", "COG0100", "COG0102", "COG0103", "COG0184", "COG0185", "COG0186", "COG0197", "COG0200", "COG0256", "COG0522", "COG0085", "COG0086", "COG0202"] # 27 ribosomal proteins and 3 RNAP subunits
		print("Using the RNAP and ribosomal marker set")
	else:
		cog_set = ["COG0012", "COG0016", "COG0018", "COG0048", "COG0049", "COG0052", "COG0080", "COG0081", "COG0085", "COG0086", "COG0087", "COG0088", "COG0090", "COG0091", "COG0092", "COG0093", "COG0094", "COG0096", "COG0097", "COG0098", "COG0099", "COG0100", "COG0102", "COG0103", "COG0124", "COG0172", "COG0184", "COG0185", "COG0186", "COG0197", "COG0200", "COG0201", "COG0202", "COG0215", "COG0256", "COG0495", "COG0522", "COG0525", "COG0533", "COG0541", "COG0552"] # all 40 proteins
		print("Using the full 40-protein marker set")

	merged_proteins = open(project+".proteins.faa", "w")
	merged = open(project+".summary.tsv", "w")
	merged.write("new_protein_name\toriginal_protein_name\tgenome_name\thit\thit_type\tprotein_length\thmm_score\thmm_alignlength\tnum_proteins_merged\tproteins_merged\tmerged_proteins_hit_coords\n")
	cog_out = open(project+".cogs.txt", "w")

	hmm_launcher(input)
	hmm_parser(input, ".domout", "all_hmmout.txt")

	for i in os.listdir(input):
		if i.endswith(".faa"):

			protein_file = os.path.join(input, i)
			gff_file = re.sub(".faa", ".gff", protein_file)
			domout = re.sub(".faa", ".domout", protein_file)
			parsed = re.sub(".faa", ".domout.parsed", protein_file)
			acc = re.sub(".faa", "", i)		

			# get a dictionary of protein sequences
			seq_handle = open(protein_file, "r")
			seq_dict = SeqIO.to_dict(SeqIO.parse(seq_handle, "fasta"))
			#orf_set = [record.id for record in seq_dict.values()]

			# parse domout file and get protein hits and coordinates
			done = {}
			for cog in cog_set:

			#	print(contig_name, orf_set)

				#print cog
				#protein2dups = defaultdict(lambda:"single_besthit")

				rnap, protein2cog, protein2acc, protein2score, protein2length, protein2category, protein2coords, protein2align_length, protein2dups = parse_domout(parsed, acc, seq_dict, cog)

				if rnap == "NAN":
					pass
				else:

					contig_name = re.sub("_\d+$", "", rnap)
					orf_set = [record.id for record in seq_dict.values() if re.sub("_\d+$", "", record.id) == contig_name]

					already_done = []
					num_proteins = defaultdict(lambda:int(1))
					prot2protlist = defaultdict(list)
					prot2loc = defaultdict(list)

					prot2protlist[rnap].append(rnap)
					id_hit1 = rnap +"|"+ cog
					range1 = protein2coords[id_hit1]

					r1 = range(range1[0], range1[1])
					meanloc1 = np.mean(range1)
					prot2loc[rnap].append(meanloc1)

					orf_set.remove(rnap)

					#print rnap
					if cog in ["COG0085", "COG0086"]:
						for d in orf_set:
							if protein2cog[d] == cog:
								id_hit2 = d +"|"+ cog
								range2 = protein2coords[id_hit2]
								r2 = range(range2[0], range2[1])
								meanloc2 = np.mean(range2)
								
								set1 = set(r1)
								inter = set1.intersection(r2)

								if int(len(inter)) > 50:
									protein2dups[id_hit2] = "false_hit"
								else:
									protein2dups[id_hit1] = "main_hit"
									protein2dups[id_hit2] = "secondary_hit"

									minrange = min(range1 + range2)
									maxrange = max(range1 + range2)
									protein2coords[id_hit1] = [minrange, maxrange]
									
									protein2align_length[id_hit1] = abs(maxrange - minrange)
									protein2length[rnap] = int(protein2length[rnap]) + int(protein2length[d])
									protein2score[rnap] = float(protein2score[rnap]) + float(protein2score[d])
									prot2protlist[rnap].append(d)
									prot2loc[rnap].append(meanloc2)
									num_proteins[id_hit1] +=1

					else:
						protein2dups[rnap +"|"+ cog]


					for item in protein2dups:
						if protein2dups[item] == "main_hit" or protein2dups[item] == "single_besthit": #or protein2dups[item] == "NEXT" or protein2dups[item] == "SECO":
							items = item.split("|")
							protein = items[0]
							hit = items[1]

							protlist = prot2protlist[protein]
							loc_list = [float(loc) for loc in prot2loc[protein]]
							index_list = [i[0] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]
							sorted_loc_list = [i[1] for i in sorted(enumerate(loc_list), key=lambda x:x[1])]

							sorted_prot_list = [protlist[index] for index in index_list]
							prot_str = ";".join(sorted_prot_list)

							#loc_str = ";".join(sorted_loc_list)
							loc_str = ";".join([str(n) for n in sorted_loc_list])
							acc = protein2acc[protein]
							final_name = re.sub("_", ".", acc) +"_"+ hit

							if hit in cog_set:
								merged.write(final_name +"\t"+ protein +"\t"+ acc +"\t"+ hit +"\t"+ protein2dups[item] +"\t"+ str(protein2length[protein]) +"\t"+ str(protein2score[protein]) +"\t"+ str(protein2align_length[item]) +"\t"+ str(num_proteins[item]) +"\t"+ prot_str +"\t"+ loc_str +"\n")

								if len(sorted_prot_list) > 1:
									newrecord = SeqRecord(Seq("", IUPAC.protein), id=final_name, name=protein+" JOINED", description=protein2acc[protein] +" JOINED PROTEIN")
									print "Joining the following proteins: " +" ".join(sorted_prot_list) +" that both have hits to: "+ cog   #newrecord.name
									for fragment in sorted_prot_list:
										subrecord = seq_dict[fragment]
										subseq = subrecord.seq
										subseq = re.sub("\*", "", str(subseq))
										newrecord.seq = newrecord.seq +""+ subseq

									final_proteins.append(newrecord)

								else:
									record = seq_dict[protein]
									record.id = final_name
									final_proteins.append(record)

	names = [i.id for i in final_proteins]
	for cog in cog_set:
		name_set = [i for i in names if cog in i]
		name_str = "\t".join(name_set)
		cog_out.write(name_str +"\n")

	final_records = []
	for seqrecord in final_proteins:
		seq = seqrecord.seq
		newseq = Seq("".join([n for n in seq if n != "*"]), IUPAC.protein)
		#print(newseq)
		newrecord = SeqRecord(newseq, id=seqrecord.id, name=seqrecord.name, description=seqrecord.description)
		final_records.append(newrecord)

	SeqIO.write(final_records, merged_proteins, "fasta")


########################################################################
##### use argparse to run through the command line options given #######
########################################################################
def main(argv=None):

	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="markerfinder- a simple script for predicting phylogenetic marker genes and formatting them for concatenated alignments \nFrank O. Aylward, Assistant Professor, Virginia Tech Department of Biological Sciences <faylward at vt dot edu>", epilog='*******************************************************************\nIf you use this tool in a publication please do not forget to cite:\nHMMER3 (DOI 10.1371/journal.pcbi.1002195)\n*******************************************************************')
	args_parser.add_argument('-i', '--input', required=True, help='Input folder of protein files (.faa extension)')
	args_parser.add_argument('-p', '--project', required=True, help='project name prefix for output files')
	args_parser.add_argument('-db', '--database', required=False, default="all", help='HMM database to use. Options are "all", "ribo", "rnap" or "ribo_rnap". See README for details')
	args_parser.add_argument('-t', '--cpus', required=False, default=str(1), help='number of cpus to use for the HMMER3 search')
	args_parser = args_parser.parse_args()

	# set up object names for input/output/database folders
	input = args_parser.input
	project = args_parser.project
	database = args_parser.database
	cpus = args_parser.cpus

	summary_file = 1
	run_program(input, project, database, cpus)

	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)

# end






























