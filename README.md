#### FILES ####
* `conf.py` : goto this python file to change all the required parameters
* `main.py`
* `alignment.py`
* `variant_call.py`


#### BOWTIE2 alignment settings ####
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
    
#### mpileup ####
```
seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
```
* chromosome(ID), 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities
* `a dot` stands for a match to the reference base on the forward strand
* `a comma` for a match on the reverse strand
* `ACGTN` for a mismatch on the forward strand
* `acgtn` for a mismatch on the reverse strand
* A pattern `\+[0-9]+[ACGTNacgtn]+` indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Here is an example of 2bp insertions on three reads:
    
    `seq2 156 A 11  .$......+2AG.+2AG.+2AGGG    <975;:<<<<<`
* Similarly, a pattern `-[0-9]+[ACGTNacgtn]+` represents a deletion from the reference. Here is an exmaple of a 4bp deletions from the reference, supported by two reads:

    `seq3 200 A 20 ,,,,,..,.-4CACC.-4CACC....,.,,.^~. ==<<<<<<<<<<<::<;2<<`
* Also at the read base column, a symbol `^` marks the start of a read segment which is a contiguous subsequence on the read separated by `N/S/H` CIGAR operations. 
* The ASCII of the character following `^` minus 33 gives the mapping quality. 
* A symbol `$` marks the end of a read segment. Start and end markers of a read are largely inspired by Phil Green's CALF format. These markers make it possible to reconstruct the read sequences from pileup.

#### VCF files ####
 * **INFO**
     * `DP`: raw read depth
     * `I16`: Auxiliary tag
     * `QS`: Auxiliary tag
     * `MQSB`: Mann-Whitney U test of mapping quality vs strand bias
     * `MQ0F`: fractrion of MQ0 reads
     
     
     