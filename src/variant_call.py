from conf import *

class VariantCall(object):
    def __init__(self, reference, bam_file, setting="DEFAULT"):
        self._reference = reference
        self._setting = setting
        self._bam_file = bam_file

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

    def _main(self):

        self._call_variants(self._bam_file)

# if __name__ == "__main__":
# 	variant_caller = VariantCall(all_reference+".fasta")
# 	variant_caller._main()