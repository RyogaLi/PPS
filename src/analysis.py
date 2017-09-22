from conf import *

def get_alignment_rate(directory):
	"""
	get alignment rate from .log file in directory
	:param directory: 
	:return: 
	"""
	percentages = []
	file_list = os.listdir(directory)
	for dirpath, dirnames, filenames in os.walk(directory):
		for filename in filenames:
			if filename.startswith("bowtie_"):
				line = subprocess.check_output(['tail', '-1', os.path.join(dirpath, filename)])
				per = float(line.split("%")[0])
				percentages.append(per)
	return percentages

def gene_count_plot(n, gc, fc):
	"""
	make a bar chart
	:param n: number of classes 
	:param y: value
	:return: 
	"""
	fig, ax = plt.subplots()
	index = np.arange(n)
	bar_width = 0.2

	opacity = 0.4
	rects1 = plt.bar(index, gc, bar_width, alpha=opacity, color='b',label='Gene count')

	rects2 = plt.bar(index + bar_width, fc, bar_width,alpha=opacity,color='r',label='Full covered gene count')
	plt.xlabel('Well')
	plt.ylabel('Count')
	plt.title('Number of genes in each well')
	plt.xticks([])
	plt.legend()

	plt.tight_layout()
	plt.savefig("./gene_count.png")

def get_full_cover(file):
	"""
	Get a dictionary of gene names which are fully covered(aligned) in vcf file
	:return: dictionary with keys = gene names; value = gene length
	"""
	with open(file, "r") as raw:
		gene_dict = {}
		ref_dict = {}
		for line in raw:
			id_line = re.search("<ID=(.+?),length=(.+?)>", line)
			if id_line:
				ref_dict[id_line.group(1)] = int(id_line.group(2))
			if "#" not in line:
				line = line.split()
				# print line[7]
				if line[0] not in gene_dict.keys():
					# grep read depth information from INFO section
					rd = re.search("DP=([0-9]+)", line[7])
					rd = rd.group(1)
					gene_dict[line[0]] = [1, int(rd)]
				else:
					# grep read depth information from INFO section
					rd = re.search("DP=([0-9]+)", line[7])
					rd = rd.group(1)
					gene_dict[line[0]][0]+=1
					gene_dict[line[0]][1]+=int(rd)
		remove_genes = gene_dict
		for key in remove_genes.keys():
			if remove_genes[key] < int(ref_dict[key]):
				del remove_genes[key]
			else:
				avg_rd = remove_genes[key][1]/ remove_genes[key][0]
				remove_genes[key][1] = avg_rd

	return remove_genes, len(gene_dict.keys()), gene_dict, ref_dict


def filter_vcf(file, gene_names):
	"""
	Filter vcf file with only genes in the gene_dictionary 
	:param gene_names: 
	:return: 
	"""
	snp_count = {}
	indel_count = {}
	read_depth = {}
	file_basename = os.path.basename(file).split(".")[0]
	# create reader and writer
	with open(file, "r") as raw_vcf:
		with open(file_basename+"_filtered.vcf", "w") as filtered:
			for line in raw_vcf:
				# eliminate header
				if line.startswith("##"):
					continue
				# remove unwanted genes
				line = line.split()
				if line[0] not in gene_names.keys() or "PDONR" in line[0]:
					continue
				read_depth[line[0]] = gene_names[line[0]][1]
				# count SNP and INDEL for each gene
				if "INDEL" in line[-3]:
					if line[0] in indel_count.keys():
						indel_count[line[0]] += 1
					else:
						indel_count[line[0]] = 1
				elif "<*>" not in line[4]:
					if line[0] in snp_count.keys():
						snp_count[line[0]] += 1
					else:
						snp_count[line[0]] = 1
				# write record to file
				filtered.write("\t".join(line)+"\n")
	return snp_count, indel_count, read_depth

#
# if __name__ == "__main__":
# 	old_ref_rate = get_alignment_rate("/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_old_ref/")
# 	combined_ref_rate = get_alignment_rate("/Users/roujia/Documents/02_dev/02_pooled_plasmid/03_PPS_dev/output_combined_ref/")
#
# 	plt.plot(range(len(old_ref_rate)), old_ref_rate, '.')
# 	plt.plot(range(len(combined_ref_rate)), combined_ref_rate, '.')
# 	plt.title("Compare alignment rate")
# 	plt.xlabel("plate")
# 	plt.ylabel("% aligned")
# 	plt.savefig("alignment_rate.png")