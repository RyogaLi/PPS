# create .sh file that can be used to submit jobs to sge
from conf import *
import conf

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Submit to SGE')
    parser.add_argument('-dir', type=str, required=True, help='input fastq path')
    args = parser.parse_args()

    input_dir = args.dir

    # take the fastq path and divide all files into groups of max_queue
    # here we run one file each time (no slave)
    list_files = os.listdir(input_dir)

    if conf.PAIRED == True:
        r1_files = []
        r2_files = []
        for file in list_files:
            if file.endswith(".fastq"):
                # fetch R1 in file name and add to r1
                if "R1" in file:
                    r1_files.append(file)
                elif "R2" in file:
                    r2_files.append(file)
        track = 0
        for r1 in r1_files:
            if track <= max_queue:
                identifier = os.path.basename(r1).split("R1")[0]
                r2 = [i for i in r2_files if identifier in i][0]
                conf.r1 = r1
                conf.r2 = r2
                # todo specify output dir and output log and error log
                cmd = "qsub -N " + identifier + " submission.sh 1> "+output
                os.system(cmd)
                track += 1
            else:
                # todo wait for all the jobs finish
                track = 0
    else:
        track = 0
        for file in list_files:
            if track <= conf.max_queue and file.endswith(".fastq"):
                conf.fastq = file
                # todo specify output dir and output log and error log
                time_stamp = strftime("%Y_%m_%d_%H_%M_%S", gmtime())
                dir_name = file.split(".")[0] + "_" + time_stamp
                output_dir = os.path.join(conf.output, dir_name)
                os.mkdir(output_dir)

                cmd = "qsub -N "+ file + " submission.sh"
                # os.system(cmd)





