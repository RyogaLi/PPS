from conf import *

class VariantCall(object):

	def __init__(self, reference, setting="DEFAULT"):
		self._reference = reference
		self._setting = setting

	def _call_variants(self, bam):
		basename =  os.path.basename(bam).split("_")[0]
		if self._setting == "DEFAULT":
			command = "samtools mpileup -uf " + self._reference + " " + bam + " | bcftools view - > " + basename + ".raw.vcf"
			# command = "samtools mpileup --ff=1024 -A -Vf " + self._reference + " " + bam + " > " + basename + ".txt"
		else:
			command = "ERROR: please provide correct setting"
			sys.exit(command)
		os.system(command)
		# os.system("bcftools view " + basename + ".raw.vcf | vcfutils.pl varFilter -D100 > " + basename + ".filtered.vcf")

	def _main(self):
		# goto each folder in output dir
		# run this inside the dir
		dir_list = os.listdir(output)
		for dir in dir_list:
			if not os.path.isdir(output+"/"+dir): continue
			os.chdir(output+dir)

			for file in os.listdir("."):
				if "_sorted.bam" in file:
					self._call_variants(file)


if __name__ == "__main__":
	variant_caller = VariantCall(reference+".fasta")
	variant_caller._main()

	# task = "bowtie2 --local -a -x %s -1 %s -2 %s -S %s; samtools view -b -o %s %s; samtools sort %s -o %s -T .temp; samtools index %s; samtools mpileup --ff=1024 -A -Vf %s %s > %s; rm %s; rm %s" % (
	# 	ref_idx, read_1, read_2, output_name + ".sam", output_name + ".bam", output_name + ".sam", output_name +  ".bam",
	# 	output_name + "_sorted.bam", output_name + "_sorted.bam", ref_file, output_name + "_sorted.bam",
	# 	output_name + ".txt", output_name + ".sam", output_name + ".bam")