#!/usr/bin/env python3.7

"""
Alignment plasmid pool sequencing analysis
"""
import os
import subprocess
import shutil

class Alignment(object):

    def __init__(self, reference, fastq, output, sh_file, log, setting="DEFAULT"):
        """
        Inicialize alignment object
        :param reference: reference files
        :param fastq: fastq file (R1 and R2 merged)
        :param output: output dir
        :param sh_file: sh file for submitting jobs
        :param log: log object
        :param setting: alignment setting
        """

        self._sample = fastq
        self._reference = reference
        # self._sub_reference = subset_reference
        self._setting = setting
        self._output = output
        self._sh_file = sh_file
        self._log = log

    def _align(self, at):
        """
        Write alignment command to sh file
        R1 and R2 are merged
        :return:
        """
        self._basename = os.path.basename(self._sample).split(".")[0]

        sam_file = os.path.join(self._output, os.path.basename(self._sample).replace(".fastq.gz", ".sam"))
        # r2_sam_file = os.path.join(output_path, os.path.basename(r2).replace(".fastq.gz", ".sam"))
        log_f = os.path.join(self._output, os.path.basename(self._sh_file).replace(".sh", ""))

        if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
            r1_cmd = f"bowtie2 -a -p 16 --local -x {self._reference} -U {self._sample} -S {sam_file}"
            # r2_cmd = f"bowtie2 -a -p 16 --local -x {self._reference} -U {r2} -S {r2_sam_file}"

        elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
            r1_cmd = f"bowtie2 -a -p 16 --local --very-sensitive-local -x {self._reference} -U {self._sample} -S {sam_file}"
            # r2_cmd = f"bowtie2 -a -p 16 --local --very-sensitive-local -x {self._reference} -U {r2} -S {r2_sam_file}"

        else:
            command = "ERROR: please provide correct setting (DEFAULT/SENSITIVE)"
            raise ValueError(command)
        time_request = f"0{at}:00:00"
        header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --job-name={self._basename}\n#SBATCH " \
                 f"--cpus-per-task=16\n#SBATCH --error={log_f}-%j.log\n#SBATCH --output={log_f}-%j.log\n"
        # write header to sh file
        with open(self._sh_file, "w") as sh:
            sh.write(header)
            sh.write(r1_cmd+"\n")
            bam_file = sam_file.replace(".sam", ".bam")
            # r2_bam_file = r2_sam_file.replace(".sam", ".bam")

            # convert sam file to a sorted bam file out put from samtools are save in corresponding log files, sterr
            sh.write(f"samtools view -bS {sam_file} > {bam_file}\n")
            sh.write(f"samtools sort {bam_file} -o {bam_file.replace('.bam', '_sorted.bam')}\n")
            # creating a bam index file
            sh.write(f"samtools index {bam_file.replace('.bam', '_sorted.bam')} "
                     f"{bam_file.replace('.bam', '_sorted.bai')}\n")
            # create vcf alignment output
            # first pileup the reads with bcftools mpileup
            sh.write(f"bcftools mpileup -f {self._reference}.fasta {bam_file.replace('.bam', '_sorted.bam')} >"
                     f" {bam_file.replace('.bam', '_raw.bcf')}\n")
            # then convert to vcf files
            sh.write(f"bcftools view -u {bam_file.replace('.bam', '_raw.bcf')} > {bam_file.replace('.bam', '_raw.vcf')}")

            # # convert sam file to a sorted bam file out put from samtools are save in corresponding log files, sterr
            # sh.write(f"samtools view -bS {r2_sam_file} > {r2_bam_file}\n")
            # sh.write(f"samtools sort {r2_bam_file} -o {r2_bam_file.replace('.bam', '_sorted.bam')}\n")
            # # creating a bam index file
            # sh.write(f"samtools index {r2_bam_file.replace('.bam', '_sorted.bam')} "
            #          f"{r2_bam_file.replace('.bam', '_sorted.bai')}\n")

    def _main(self, at):
        """
        Write sh file and submit sh file for alignment
        :param r1:
        :param r2:
        :return:
        """
        self._align(at)
        os.system(f"chmod 755 {self._sh_file}")
        # submit this to the cluster
        sub_cmd = ["sbatch", str(self._sh_file)]
        self._log.debug(sub_cmd)
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
        # log sample name and job id
        self._log.info(f"Sample {self._basename}: job id - {job_id}")
        return job_id.split()[-1]

