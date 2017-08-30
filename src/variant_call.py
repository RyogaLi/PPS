from conf import *

class VariantCall(object):


	def __init__(self, reference, setting="DEFAULT"):
		self._reference = reference
		self._setting = setting

	def _call_variants(self, bam):
		basename =  os.path.basename(bam).split(".")[0]
		if self._setting == "DEFAULT":
			command = "samtools mpileup -uf " + self._reference + " " + bam + " | bcftools view - > " + basename + ".raw.vcf"
			self._raw_vcf = basename + ".raw.vcf"
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

	def _filter_vcf(self):
		# todo filter all the reads with 1 rd
		# todo filter all the genes that are not fully covered
		# todo separate into snp and indel

		with open(self._raw_vcf, "r") as raw:
			gene_dict = {}
			ref_dict = {}
			for line in raw:
				id_line = re.search("<ID=(.+?),length=(.+?)>", line)
				if id_line:
					ref_dict[id_line.group(1)] = int(id_line.group(2))
				if "#" not in line:
					if "INDEL" in line: continue
					line = line.strip().split("\t")
					# INFO = line[7].split(";")
					# print INFO
					# dp = int(INFO[0].split("=")[1])
					if line[0] not in gene_dict.keys():
						# if dp == 0: continue
						gene_dict[line[0]] = 1
					else:
						gene_dict[line[0]] +=1

			for key in gene_dict.keys():

				if gene_dict[key] != int(ref_dict[key]):
					print key
					print gene_dict[key]
					print ref_dict[key]
					del gene_dict[key]

		return gene_dict


	def _main(self):
		# goto each folder in output dir
		# run this inside the dir
		read_depth = {}
		dir_list = os.listdir(output)
		for dir in dir_list:
			if not os.path.isdir(output+"/"+dir): continue
			os.chdir(output+dir)

			for file in os.listdir("."):
				if "_sorted.bam" in file:
					# call variant
					self._call_variants(file)

					# summarize read depth for all genes
					#read_depth = self._analyze_rd()
					print self._filter_vcf()
					break


if __name__ == "__main__":
	variant_caller = VariantCall(reference+".fasta")
	variant_caller._main()
	# test parse

	# task = "bowtie2 --local -a -x %s -1 %s -2 %s -S %s; samtools view -b -o %s %s; samtools sort %s -o %s -T .temp; samtools index %s; samtools mpileup --ff=1024 -A -Vf %s %s > %s; rm %s; rm %s" % (
	# 	ref_idx, read_1, read_2, output_name + ".sam", output_name + ".bam", output_name + ".sam", output_name +  ".bam",
	# 	output_name + "_sorted.bam", output_name + "_sorted.bam", ref_file, output_name + "_sorted.bam",
	# 	output_name + ".txt", output_name + ".sam", output_name + ".bam")