from conf import *

class VariantCall(object):
	def __init__(self, reference, setting="DEFAULT"):
		self._reference = reference
		self._setting = setting


	def _call_variants(self, bam):
		self._basename = os.path.basename(bam).split(".")[0]
		if self._setting == "DEFAULT":
			command = "samtools mpileup -uf " + self._reference + " " + bam + " | bcftools view - > " + self._basename + ".raw.vcf"
			self._raw_vcf = self._basename + ".raw.vcf"
			# command = "samtools mpileup --ff=1024 -A -Vf " + self._reference + " " + bam + " > " + basename + ".txt"
		else:
			command = "ERROR: please provide correct setting"
			sys.exit(command)
		os.system(command)
		# os.system("bcftools view " + basename + ".raw.vcf | vcfutils.pl varFilter -D100 > " + basename + ".filtered.vcf")
		return self._raw_vcf

	def _analyze_rd(self):
		gene_list = {}
		vcf_reader = vcf.Reader(open(self._raw_vcf))
		for record in vcf_reader:
			depth = int(record.INFO["DP"])
			gene = record.CHROM
			if gene not in gene_list.keys():
				gene_list[gene]=[depth]
			else:
				gene_list[gene].append(depth)
		total = 0
		count = 0
		for key in gene_list.keys():
			total += 1
			gene_list[key]= sum(gene_list[key])/len(gene_list[key])
			if gene_list[key] < 1 or gene_list[key] >1000: # set threshold
				count += 1
				print key
				print gene_list[key]
				del gene_list[key]

		# save info to file
		with open("read_depth.txt", "w") as rd_file:
			rd_file.write("gene_name\tread_depth\n")
			for key in gene_list.keys():
				rd_file.write(key + "\t" + str(gene_list[key]) + "\n")
		# plot read depth
		rd = gene_list.values()
		rd.sort()
		plt.plot(range(len(rd)), rd, '.')
		plt.xlabel("Genes")
		plt.ylabel("read depth")
		plt.grid(True)
		plt.savefig("./Read_depth.png")
		plt.close()

		return gene_list

	def _get_full_cover(self):
		"""
		Get a dictionary of gene names which are fully covered(aligned) in vcf file
		:return: dictionary with keys = gene names value = gene length
		"""
		with open(self._raw_vcf, "r") as raw:
			gene_dict = {}
			ref_dict = {}
			total_gene_count= []
			for line in raw:
				id_line = re.search("<ID=(.+?),length=(.+?)>", line)
				if id_line:
					ref_dict[id_line.group(1)] = int(id_line.group(2))
				if "#" not in line:
					# if "INDEL" in line: continue
					line = line.split()
					if line[0] not in gene_dict.keys():
						gene_dict[line[0]] = 1
					else:
						gene_dict[line[0]] +=1

					if line[0] not in total_gene_count:
						total_gene_count.append(line[0])

			for key in gene_dict.keys():
				if gene_dict[key] < int(ref_dict[key]):
					del gene_dict[key]

		return gene_dict, len(total_gene_count)

	def _filter_vcf(self, gene_names):
		"""
		Filter vcf file with only genes in the gene_dictionary 
		:param gene_names: 
		:return: 
		"""
		snp_count = {}
		indel_count = {}
		# create reader and writer
		with open(self._raw_vcf, "r") as raw_vcf:
			with open(self._basename+"_filtered.vcf", "w") as filtered:
				for line in raw_vcf:
					# eliminate header
					if line.startswith("##"):
						continue
					# remove unwanted genes
					line = line.split()
					if line[0] not in gene_names.keys():
						continue
					# count SNP and INDEL for each gene

					if "INDEL" in line[-3]:
						if line[0] in indel_count.keys():
							indel_count[line[0]] += 1
						else:
							indel_count[line[0]] = 1
					elif line[4] != "<*>":
						if line[0] in snp_count.keys():
							snp_count[line[0]] += 1
						else:
							snp_count[line[0]] = 1

					# write record to file
					filtered.write("\t".join(line)+"\n")
		return snp_count, indel_count

	def _gene_count_plot(self, n, gc, fc):
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
		plt.xticks(index + bar_width / 2)
		plt.legend()

		plt.tight_layout()
		plt.show()

	def _main(self):
		# goto each folder in output dir
		# run this inside the dir
		read_depth = {}
		gc = []
		fc = []
		total_files = 0
		dir_list = os.listdir(output)
		for dir in dir_list:
			if not os.path.isdir(output+"/"+dir): continue
			os.chdir(output+dir)
			for file in os.listdir("."):
				if "_sorted.bam" in file:
					total_files += 1
					# call variant
					self._call_variants(file)
					# get genes that are fully covered by alignment
					full_cover, total = self._get_full_cover()
					snp, indel = self._filter_vcf(full_cover)
					gc.append(total)
					fc.append(len(full_cover.keys()))
					if total_files >= 3: break
		self._gene_count_plot(total_files, gc, fc)

#
if __name__ == "__main__":
	variant_caller = VariantCall(reference+".fasta")
	variant_caller._main()
	# test plot
	# a = [2,3,4]
	# b = [3,4,5]
	# n = 3
	# variant_caller._gene_count_plot(n, a, b)

	# task = "bowtie2 --local -a -x %s -1 %s -2 %s -S %s; samtools view -b -o %s %s; samtools sort %s -o %s -T .temp; samtools index %s; samtools mpileup --ff=1024 -A -Vf %s %s > %s; rm %s; rm %s" % (
	# 	ref_idx, read_1, read_2, output_name + ".sam", output_name + ".bam", output_name + ".sam", output_name +  ".bam",
	# 	output_name + "_sorted.bam", output_name + "_sorted.bam", ref_file, output_name + "_sorted.bam",
	# 	output_name + ".txt", output_name + ".sam", output_name + ".bam")