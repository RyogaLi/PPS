import numpy as np
import os

def process_seq_file(seq_file):
	"""
	remove all the dashes in the gene name
	remove SUP sequence if the same gene is in HIP
	save the duplicated SUP sequence to another file
	:param seq_file: all_sequence.txt
	:return: 
	"""
	with open(seq_file, "r") as seq, open("02_all_sequence_HIP.txt", "w") as out, open("03_all_sequence_other.txt", "w") as rep_file:
		out.write("orf_name\tsource\tcds_seq\taa_seq\n")
		rep_file.write("orf_name\tsource\tcds_seq\taa_seq\n")
		for line in seq:
			if "orf_name" in line: continue # skip colname
			line = line.strip().split("\t")
			ORF_name = line[0]
			source = line[1]
			if "-" in ORF_name: # remove all the dashes in orf_name
				ORF_name = "".join(ORF_name.split("-"))
			if source == "HIP":
				out.write(ORF_name+"\t"+source+"\t"+line[2]+"\t"+line[3]+"\n")
			else:
				rep_file.write(ORF_name+"\t"+source+"\t"+line[2]+"\t"+line[3]+"\n")


def read_file_to_matrix(filename):
	output = []
	with open(filename, "r") as file_in:
		for line in file_in:
			line = line.strip().split()
			output.append(line)
	return np.asmatrix(output)

def make_fasta(input_file):
	"""
	make a reference (fasta file)
	:param input_file: 
	:return: 
	"""
	with open(input_file, "r") as fp, open("ORF_reference.fasta", "w") as ref:
		for line in fp:
			if "source" in line: continue
			line = line.strip().split()
			ref.write(">"+line[0]+"_"+line[1]+"\n")
			ref.write(line[2]+"\n")


if __name__ == "__main__":
	# seq_file = "/Users/roujia/Documents/02_dev/02_pooled_plasmid/02_data/all_sequence.txt"
	# process_seq_file(seq_file)

	# hip = read_file_to_matrix("./02_all_sequence_HIP.txt")
	# other = read_file_to_matrix("./03_all_sequence_other.txt")
	# fasta = "ORF_ref.fasta"
	# os.system('cat ./02_all_sequence_HIP.txt > ./04_removed_dup.txt')
	# command = "awk 'FNR==NR{a[$1];next}!($1 in a)' ./02_all_sequence_HIP.txt ./03_all_sequence_other.txt >> 04_removed_dup.txt"
	# os.system(command)
	make_fasta("./04_removed_dup.txt")