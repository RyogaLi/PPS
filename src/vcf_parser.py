import vcf

class VcfFileParser(object):

    def __init__(self, vcf):
        self._vcf = vcf

    def _parse(self):
        gene_variant_counter = {}
        vcf_reader = vcf.Reader(self._vcf)
        indel_file = self._vcf.split(".")[0]+"_indel.vcf"
        indel_writer = vcf.Writer(indel_file)
        for record in vcf_reader:
            if record.ALT == "<*>": continue
            if "INDEL" in record.INFO.keys():
                indel_writer.write_record(record)
