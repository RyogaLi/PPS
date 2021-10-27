#!/usr/bin/env python3.7

"""
Alignment plasmid pool sequencing analysis
"""
import os
import re
import subprocess
import shutil

class Alignment(object):

    def __init__(self, ref_dir, mode, fastq, output, sh_file, log, setting="DEFAULT"):
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
        # self._reference = reference
        self._setting = setting
        self._output = output
        self._sh_file = sh_file
        self._log = log
        self._basename = os.path.basename(self._sample).split("_")[0]
        self._mode = mode
        self._ref_dir = ref_dir

    def _align(self, reference, suffix):
        """
        Write alignment command to sh file
        R1 and R2 are merged
        :return:
        """
        sam_file = os.path.join(self._output, os.path.basename(self._sample).replace(".fastq.gz", f"{suffix}.sam"))

        if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
            r1_cmd = f"bowtie2 -a -p 4 --local -x {reference} -U {self._sample} -S {sam_file}"
            # r2_cmd = f"bowtie2 -a -p 16 --local -x {self._reference} -U {r2} -S {r2_sam_file}"

        elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
            r1_cmd = f"bowtie2 -a -p 4 --local --very-sensitive-local -x {reference} -U {self._sample} -S {sam_file}"
            # r2_cmd = f"bowtie2 -a -p 16 --local --very-sensitive-local -x {self._reference} -U {r2} -S {r2_sam_file}"

        else:
            command = "ERROR: please provide correct setting (DEFAULT/SENSITIVE)"
            raise ValueError(command)

        # write header to sh file
        with open(self._sh_file, "a") as sh:
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
            sh.write(f"bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,"
                     f"INFO/ADF,INFO/ADR -f {reference}.fasta {bam_file.replace('.bam', '_sorted.bam')} >"
                     f" {bam_file.replace('.bam', '_raw.bcf')}\n")
            # then convert to vcf files
            sh.write(f"bcftools view -u {bam_file.replace('.bam', '_raw.bcf')} > {bam_file.replace('.bam', '_raw.vcf')}\n")
            
            # get vcf file with variants only
            sh.write(f"bcftools call -mAv --ploidy 1 {bam_file.replace('.bam', '_raw.bcf')} >"
                     f" {bam_file.replace('.bam', '_variants.vcf')}\n\n")
            # # convert sam file to a sorted bam file out put from samtools are save in corresponding log files, sterr
            # sh.write(f"samtools view -bS {r2_sam_file} > {r2_bam_file}\n")
            # sh.write(f"samtools sort {r2_bam_file} -o {r2_bam_file.replace('.bam', '_sorted.bam')}\n")
            # # creating a bam index file
            # sh.write(f"samtools index {r2_bam_file.replace('.bam', '_sorted.bam')} "
            #          f"{r2_bam_file.replace('.bam', '_sorted.bai')}\n")

    def main(self, at):
        """
        Write sh file and submit sh file for alignment
        :param at: alignment time when submiting jobs to the cluster
        :return:
        """
        log_f = os.path.join(self._output, os.path.basename(self._sh_file).replace(".sh", ""))
        time_request = f"0{at}:00:00"
        header = f"#!/bin/bash\n#SBATCH --time={time_request}\n#SBATCH --job-name={self._basename}\n#SBATCH " \
                 f"--cpus-per-task=4\n#SBATCH --error={log_f}-%j.log\n#SBATCH --output={log_f}-%j.log\n"
        # write header to sh file
        with open(self._sh_file, "w") as sh:
            sh.write(header)

        if self._mode == "yeast":
            # get reference for yeast
            # for each plate, align to all_orfs, plate_orfs, subset_orfs
            # files labeld with HIP means they are in the HIP group
            # files labeled with Sup01-Sup03 means they are in the supp-PROTGEN group
            # files labeled with Sup04+ means they are in the supp-SGD group
            # all orfs
            all_orfs = os.path.join(self._ref_dir, "all_seq")
            plate_orfs = os.path.join(self._ref_dir, self._basename)
            regexp = re.compile(r"Sup0[1-3]")
            if "HIP" in self._basename:
                sub_set_orfs = os.path.join(self._ref_dir, "hip_all")
            elif regexp.search(self._basename):  # small orfs
                sub_set_orfs = os.path.join(self._ref_dir, "PROTGEN_all")
            else:
                sub_set_orfs = os.path.join(self._ref_dir, "SGD_all")

            # in the final version, we only align orfs to the targeted plate ORFs
            # self._align(all_orfs_backbone, "_allwithbackbone")
            # self._align(all_orfs, "_allORFs")
            # self._align(sub_set_orfs, "_subsetORFs")
            self._align(plate_orfs, "_plateORFs")

        else:  # human samples
            # get group name for this sample
            match = re.search(".+(G[0-9]+)", self._basename)
            if match:
                group_name = match.group(1)
            else:
                self._log.info(self._basename)
                raise ValueError("no group ID specified")

            group_spec_orfs = os.path.join(self._ref_dir, f"group_ref_{group_name}")
            #self._align(all_orfs, "_all_orfs")
            self._align(group_spec_orfs, "_group_spec_orfs")

        os.system(f"chmod 755 {self._sh_file}")
        # submit this to the cluster
        sub_cmd = ["sbatch", str(self._sh_file)]
        self._log.debug(sub_cmd)
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
        # log sample name and job id
        self._log.info(f"Sample {self._basename}: job id - {job_id}")

        return job_id.split()[-1]

