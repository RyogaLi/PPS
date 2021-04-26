from conf import *


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
				ref_dict[id_line.group(1)] = int(id_line.group(2)) # assign gene length to geneID
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
				remove_genes[key].append(avg_rd)

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
				if "INDEL" in line[-3]: # an INDEL
					if line[0] in indel_count.keys():
						indel_count[line[0]] += 1
					else:
						indel_count[line[0]] = 1
				# elif "<*>" not in line[4]: # SNP
				elif "<*>" not in line[4]:

					# for each SNP, find out position and ALT
					if line[0] in snp_count.keys():
						snp_count[line[0]].append((int(line[1]),line[4][0]))
					else:
						snp_count[line[0]] = [(int(line[1]),line[4][0])]
				# write record to file

				filtered.write("\t".join(line)+"\n")
	# print len(snp_count["70364"])
	return snp_count, indel_count, read_depth

def remove_synonymous(snp_dict, dna_seq):
	"""
	for each mutation in snp_dict, find out if the mutation is synomounos 
	:param snp_dict: 
	:param dna_seq: 
	:return: 
	"""
	non_syn_dict = {}
	for protein in snp_dict.keys():
		# get dna seq from dna_seq
		# print snp_dict[protein]
		dna = dna_seq[protein]
		# print protein
		# print dna
		# print dna
		altered_dna = list(dna)
		# print altered_dna

		# alter all the bp in original dna seq
		for pos in snp_dict[protein]:
			altered_dna[pos[0]-1] = pos[1]
		altered_dna = "".join(altered_dna)
		dna = Seq(dna, generic_dna)
		ref_protein = dna.translate()
		# print ref_protein
		altered_dna = Seq(altered_dna, generic_dna)
		altered_protein = altered_dna.translate()
		# print altered_protein
		non_syn = []
		# find position that are different

		for i in range(len(altered_protein)):
			if altered_protein[i] != ref_protein[i]:
				if i == 0:
					dna_range=range(i,(i+1)*3)
				else:
					dna_range = range(i*3, (i+1)*3)
				snp = dict(snp_dict[protein])
				# print "test"+str(len(snp.keys()))
				# for all the snp in that range, add them to non-syn
				for pos in dna_range:
					try:
						# print (pos, snp[pos])
						non_syn.append((pos, snp[pos]))
					except Exception as e:
						# print("exc")
						pass
		non_syn_dict[protein] = non_syn
	# print len(non_syn_dict["70364"])
	return non_syn_dict