## Evolutionary contingency in human H3N2 influenza neuraminidase   

### Dependencies    
* [PEAR](https://github.com/tseemann/PEAR)

### Input files   
* [./Fasta/SD93_mutlib_ref.fa](./Fasta/SD93_mutlib_ref.fa): Reference sequence for SD93 mutant library
* Raw read files in fastq format from XXXXXXX

### Identifying permissive mutations for N387K
1. Merge overlapping paired-end reads
```pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]```
    - Output files should be named as:
      - SD93_ipt.fastq for input library
      - SD93_rep1.fastq for post-passaged library (replicate 1)
      - SD93_rep2.fastq for post-passaged library (replicate 2)

### Identifying permissive mutations for N336H
1. Merge overlapping paired-end reads
