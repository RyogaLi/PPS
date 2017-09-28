from conf import *

def count_gene(summary, col):
	"""
	for each plate in summary, count all the genes
	:param summary: summary file
	:param col: column that contains gene name
	:return: 
	"""
	all_plates = {}
	with open(summary, "r") as total:
		for line in total:
			if "plate" in line: continue
			line = line.split(",")
			plate_name = line[0]
			if plate_name not in all_plates.keys():
				all_plates[plate_name] = [line[col]]
			else:
				all_plates[plate_name].append(line[col])
	return all_plates


def find_overlap(d1_ref, d2):
	"""
	find overlap between values with the same key
	:param d1: 
	:param d2: 
	:return: 
	"""
	summary = {}
	for key in d1_ref.keys():
		ref_genes = set(d1_ref[key])
		genes = set(d2[key])
		# number of genes we got
		recovered_genes = len(list(ref_genes&genes))
		# number of all the genes
		all_genes = len(ref_genes)

		summary[key] = [recovered_genes,all_genes]
	return summary


def reads_count(fastq_file_list):
	"""
	count number of reads in each fastq file
	:param fastq: path to fastq file
	:return: 
	"""
	rc = {}
	for line in fastq_file_list:
		cmd = "cat "+line.strip()+" | wc -l"
		output = subprocess.Popen(cmd, stdout=subprocess.PIPE)
		temp = int(output.stdout.read())
		rc[line] = temp

	return rc

if __name__ == '__main__':
	ref_summary_hip = ""
	ref_summary_sup = ""
	summary = ""

	hip_genes = count_gene(ref_summary_hip,0)
	sup_genes = count_gene(ref_summary_sup,0)
	found_genes = count_gene(summary,0)

	# concatenate hip and sup
	concatenate = {}
	# find overlap
	overlap_dict = find_overlap(concatenate, found_genes)

	fastq_files = ""

	all_reads = reads_count(fastq_files)

	# plot from (overlap and reads count)