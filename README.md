## Evolutionary contingency in human H3N2 influenza neuraminidase   

### Dependencies    
* [PEAR](https://github.com/tseemann/PEAR)

### Input files   
* [./Fasta/SD93_mutlib_ref.fa](./Fasta/SD93_mutlib_ref.fa): Reference sequence for SD93 mutant library
* Raw read files in fastq format from XXXXXXX

### Identifying permissive mutations for N387K
1. Merge overlapping paired-end reads   
```pear -f [FASTQ FILE FOR FORWARD READ] -r [FASTQ FILE FOR FORWARD READ] -o [OUTPUT FASTQ FILE]```
    - Output files should be placed in the fastq/ folder and named as:
      - SD93_ipt.assembled.fastq for input library
      - SD93_rep1.assembled.fastq for post-passaged library (replicate 1)
      - SD93_rep2.assembled.fastq for post-passaged library (replicate 2)

2. Count mutants from fastq files   
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
