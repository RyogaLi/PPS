#!/usr/bin/env python3.7

"""
Alignment plasmid pool sequencing analysis
"""
import os
import subprocess
import sys
import logging.config
import shutil

class Alignment(object):

    def __init__(self, all_reference, fastq_id, output, sh_file, log, setting="DEFAULT", paired=False):
        self._sample_id = fastq_id
        self._reference = all_reference
        # self._sub_reference = subset_reference
        self._setting = setting
        self._paird = paired
        self._output = output
        self._sh_file = sh_file
        self._log = log

    def _align(self, r1, r2, at):
        """
        Write alignment command to sh file
        :param r1: read one of the
        :param r2:
        :return:
        """
        self._basename = os.path.basename(r1).split(".")[0]
        # create a dir for each alignment
        output_path = os.path.join(self._output, self._basename)
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.makedirs(output_path)

        r1_sam_file = os.path.join(output_path, os.path.basename(r1).replace(".fastq.gz", ".sam"))
        r2_sam_file = os.path.join(output_path, os.path.basename(r2).replace(".fastq.gz", ".sam"))
        log_f = os.path.join(output_path, os.path.basename(self._sh_file).replace(".sh", ".log"))

        if self._setting == "DEFAULT": # default bowtie2 settings for alignment, more info in README
            r1_cmd = f"bowtie2 -a -p 16 --local -x {self._reference} -U {r1} -S {r1_sam_file}"
            r2_cmd = f"bowtie2 -a -p 16 --local -x {self._reference} -U {r2} -S {r2_sam_file}"

        elif self._setting == "SENSITIVE": # strict bowtie2 settings for alignment, more info in README
            r1_cmd = f"bowtie2 -a -p 16 --local --very-sensitive-local -x {self._reference} -U {r1} -S {r1_sam_file}"
            r2_cmd = f"bowtie2 -a -p 16 --local --very-sensitive-local -x {self._reference} -U {r2} -S {r2_sam_file}"

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
            sh.write(r2_cmd+"\n")
            r1_bam_file = r1_sam_file.replace(".sam", ".bam")
            r2_bam_file = r2_sam_file.replace(".sam", ".bam")

            # convert sam file to a sorted bam file out put from samtools are save in corresponding log files, sterr
            sh.write(f"samtools view -bS {r1_sam_file} > {r1_bam_file}\n")
            sh.write(f"samtools sort {r1_bam_file} -o {r1_bam_file.replace('.bam', '_sorted.bam')}\n")
            # creating a bam index file
            sh.write(f"samtools index {r1_bam_file.replace('.bam', '_sorted.bam')} "
                     f"{r1_bam_file.replace('.bam', '_sorted.bai')}\n")

            # convert sam file to a sorted bam file out put from samtools are save in corresponding log files, sterr
            sh.write(f"samtools view -bS {r2_sam_file} > {r2_bam_file}\n")
            sh.write(f"samtools sort {r2_bam_file} -o {r2_bam_file.replace('.bam', '_sorted.bam')}\n")
            # creating a bam index file
            sh.write(f"samtools index {r2_bam_file.replace('.bam', '_sorted.bam')} "
                     f"{r2_bam_file.replace('.bam', '_sorted.bai')}\n")

    def _main(self, r1, r2, at=12):
        """
        Write sh file and submit sh file for alignment
        :param r1:
        :param r2:
        :return:
        """
        self._align(r1, r2, at)
        os.system(f"chmod 755 {self._sh_file}")
        # submit this to the cluster
        sub_cmd = ["sbatch", str(self._sh_file)]
        self._log.debug(sub_cmd)
        job = subprocess.run(sub_cmd, stdout=subprocess.PIPE)
        job_id = job.stdout.decode("utf-8").strip()
        # log sample name and job id
        self._log.info(f"Sample {self._basename}: job id - {job_id}")
        return job_id.split()[-1]

