#!/usr/bin/env python3.7

# Author: Roujia Li
# email: Roujia.li@mail.utoronto.ca

import sys

sys.path.append('..')

import os
import pandas as pd
import re
import subprocess
import time


def parse_jobs_galen(job_list, logger):
    """
    Galen uses slurm scheduler, different from BC and DC
    return true if all the jobs in job list finished
    else wait for 10 mins and return how man jobs are running and queued
    job_list: list of job ids
    logger: logging object
    """
    qstat_cmd = ["squeue", "-j", ",".join(job_list)]
    job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    qstat_out = job.stdout.decode("utf-8", errors="replace")
    qstat_err = job.stderr.decode("utf-8", errors="replace")
    logger.debug(qstat_out)
    logger.debug(qstat_err)
    while True:
        running_jobs, queued_jobs, completed_jobs = [], [], []
        if qstat_out != "":
            # make df
            qstat_df = pd.DataFrame([i.split() for i in qstat_out.split("\n")])
            qstat_df = qstat_df.rename(columns=qstat_df.iloc[0])
            qstat_df = qstat_df.drop(qstat_df.index[1])
            logger.debug(qstat_df)
            # get all active job ID
            running_jobs = qstat_df[qstat_df["ST"] == "R"]["JOBID"].tolist()
            # in queue
            queued_jobs = qstat_df[qstat_df["ST"] == "PD"]["JOBID"].tolist()
            # completing
            completed_jobs = qstat_df[qstat_df["ST"] == "CG"]["JOBID"].tolist()
            logger.debug(running_jobs)
            logger.debug(queued_jobs)

        logger.info(f"{len(queued_jobs)} jobs queued")
        logger.info(f"{len(running_jobs)} jobs running")

        final_list = list(set(running_jobs + queued_jobs))
        if final_list == []:
            return True
        else:
            # check in 10min
            time.sleep(600)
            job_list = final_list
            qstat_cmd = ["squeue", "-j", ",".join(job_list)]
            job = subprocess.run(qstat_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            qstat_out = job.stdout.decode("utf-8", errors="replace")
            qstat_err = job.stderr.decode("utf-8", errors="replace")
