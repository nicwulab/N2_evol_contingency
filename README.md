## Evolutionary contingency in human H3N2 influenza neuraminidase   
This README describes the deep sequencing analysis in:   
[Prevalence and mechanisms of evolutionary contingency in human influenza H3N2 neuraminidase](https://www.biorxiv.org/content/10.1101/2022.02.24.481718v1)

### Contents
* [Identifying permissive mutations for N387K](#Identifying-permissive-mutations-for-N387K)
* [Identifying permissive mutations for N336H](#Identifying-permissive-mutations-for-N336H)

### Dependencies    
* [Python3](https://www.python.org/) 
* [PEAR](https://github.com/tseemann/PEAR)
* [BioPython](https://github.com/biopython/biopython)
* [Distance](https://pypi.org/project/Distance/)
* [R](https://www.r-project.org/)

### Input files   
* [./Fasta/SD93_mutlib_ref.fa](./Fasta/SD93_mutlib_ref.fa): Reference sequence for SD93 mutant library
* [./Fasta/Bil69_mutlib_ref.fa](./Fasta/Bil69_mutlib_ref.fa): Reference sequence for Bil69 mutant library
* Raw read files in fastq format from NIH SRA database [BioProject PRJNA790468](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA790468)

### Identifying permissive mutations for N387K
1. Merge overlapping paired-end reads using [PEAR](https://github.com/tseemann/PEAR)  
```pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]```
    - Output files should be placed in the fastq/ folder and named as:
      - SD93_ipt.assembled.fastq for input library
      - SD93_rep1.assembled.fastq for post-passaged library (replicate 1)
      - SD93_rep2.assembled.fastq for post-passaged library (replicate 2)

2. Count variants from fastq files   
```python3 script/SD93_NA_fastq2count.py```
    - Input files:
      - [./Fasta/SD93_mutlib_ref.fa](./Fasta/SD93_mutlib_ref.fa)
      - fastq/SD93_ipt.assembled.fastq
      - fastq/SD93_rep1.assembled.fastq
      - fastq/SD93_rep2.assembled.fastq
    - Output file:
      - [./result/SD93_MultiMutLib.tsv](./result/SD93_MultiMutLib.tsv)

3. Filtering low count variants   
```python3 script/SD93_NA_filter.py```
    - Input file:
      - [./result/SD93_MultiMutLib.tsv](./result/SD93_MultiMutLib.tsv)
    - Output file:
      - [./result/SD93_MultiMutLib_filtered.tsv](./result/SD93_MultiMutLib_filtered.tsv)

4. Plot post-passaged frequency of each variant   
```Rscript script/SD93_NA_plot_freq.R```
   - Input file:
     - [./result/SD93_MultiMutLib_filtered.tsv](./result/SD93_MultiMutLib_filtered.tsv)
   - Output file:
     - [./graph/SD93_mutlib_rep_compare.png](./graph/SD93_mutlib_rep_compare.png)

### Identifying permissive mutations for N336H
1. Merge overlapping paired-end reads   
```pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]```
    - Output files should be placed in the fastq/ folder and named as:
      - Bil69_ipt.assembled.fastq for input library
      - Bil69_rep1.assembled.fastq for post-passaged library (replicate 1)
      - Bil69_rep2.assembled.fastq for post-passaged library (replicate 2)

2. Compute enrichment value for each variant from fastq files   
```python3 script/Bil69_NA_fastq2enrich.py```
    - Input files:
      - [./Fasta/Bil69_mutlib_ref.fa](./Fasta/Bil69_mutlib_ref.fa)
      - fastq/Bil69_ipt.assembled.fastq
      - fastq/Bil69_rep1.assembled.fastq
      - fastq/Bil69_rep2.assembled.fastq
    - Output files:
      - [./result/Bil69_MultiMutLib.tsv](./result/Bil69_MultiMutLib.tsv)

3. Filter irrelevant variants   
```python3 script/Bil69_NA_filter.py```
    - Input file:
      - [./result/Bil69_MultiMutLib.tsv](./result/Bil69_MultiMutLib.tsv)
    - Output file:
      - [./result/Bil69_MultiMutLib_filtered.tsv](./result/Bil69_MultiMutLib_filtered.tsv)

4. Plot enrichment of each variant   
```Rscript script/Bil69_NA_plot_enrich.R```
    - Input file:
      - [./result/Bil69_MultiMutLib.tsv](./result/Bil69_MultiMutLib.tsv)
    - Output file:
      - [./graph/Bil69_mutlib_rep_compare.png](./graph/Bil69_mutlib_rep_compare.png)
