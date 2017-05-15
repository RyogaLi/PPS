from conf import *

class VariantCall(object):

	def __init__(self, reference, setting="DEFAULT"):
		self._reference = reference
		self._setting = setting

	def _call_variants(self, bam):
		basename =  os.path.basename(bam).split("_")[0]
		if self._setting == "DEFAULT":
			command = "samtools mpileup -uf " + self._reference + " " + bam + " | bcftools view - > " + basename + ".raw.bcf 2> variant_call_raw_bcf.log"
		else:
			command = "ERROR: please provide correct setting"
			sys.exit(command)
		os.system(command)
		os.system("bcftools view " + basename + ".raw.bcf | vcfutils.pl varFilter -D100 > " + basename + ".filtered.vcf")

	def _main(self):
		# goto each folder in output dir
		# run this inside the dir
		dir_list = os.listdir(output)
		for dir in dir_list:
			os.chdir(output+dir)
			for file in os.listdir("."):
				if "_sorted.bam" in file:
					self._call_variants(file)


if __name__ == "__main__":
	variant_caller = VariantCall(reference+".fasta")
	variant_caller._main()