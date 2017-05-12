#### FILES ####
* `main.py`
* `alignment.py`
* `variant_call.py`


#### BOWTIE2 alignment settings
* **DEFAULT**
    * **command**: `bowtie2 -a -x reference -U fastq -S sam_file`
        * end to end alignment
        * `-a`: no upper limit for alignment
* **SENSITIVE**
    * **command**: `bowtie2 -a -p 20 -x ORF_with_pDONR --local --very-sensitive-local -U fastq -S sam_file`
        * `-a`: no upper limit of alignment  
        * `-p 20`: 20 thread 
        * `--local`: local alignment 
        * `--very-sensitive-local`: `-D 20 -R 3 -i S,1,0,50 -N 0 -L 20`
        * `-D 20`: 
        * `-R 3`: the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds
        * `-i S,1,0.50`:  `f(read_length) = 1 + 0.5 * sqrt(read_length)` 
        * `-N 0`: zero mismaches
        * `-L 20`: Sets the length of the seed substrings to align during multiseed alignment. Smaller values make alignment slower but more sensitive. Default: the --sensitive preset is used by default, which sets -L to 20 both in --end-to-end mode and in --local mode.