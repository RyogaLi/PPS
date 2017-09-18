## Dependencies
* [python version 2.7.10](https://www.python.org/downloads/)

* [Bowtie2 version 2.3.2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/)

* [samtools version 1.4](https://sourceforge.net/projects/samtools/files/)

## Input
* Reference file: `reference.fasta`
* For paired fastq files: `*_R1.fastq` and `*_R2.fastq` 
* For unpaired fastq files: `*.fastq`
* To combine unpaired r1 and r2 fastq files, please provide path to fastq files in `conf.py` run:
    
    ```
    python src/sup.py combine
    ```
    Combined fastq files will be saved into the same directory as old fastq files
 

## Run
0. If you haven't build index for the fasta file, please run the following command to build index files
    
    ```
    bowtie2-build [options]* path_to_your_ref.fasta reference_name
    ```
    
1. Please goto `./src/conf.py` to change all the required input parameters
2. 

## Output


## Test Data


## Notes
* Please goto [wiki page](https://github.com/RyogaLi/PPS/wiki) for more information