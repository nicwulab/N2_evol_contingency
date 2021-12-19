#!/usr/bin/python
import os
import sys
import distance
from Bio import SeqIO
from collections import Counter

def translation(seq):
  dnamap = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"_", "TAG":"_",
    "TGT":"C", "TGC":"C", "TGA":"_", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
  pep = []
  i = 0
  while i < len(seq):
    codon = seq[i:i+3]
    aa = dnamap[codon]
    pep.append(aa)
    i = i + 3
  pep = ''.join(pep)
  return pep

def ProcessMultilib(merged_file, trim_5, trim_3):
  print ("Reading %s" % merged_file)
  records = SeqIO.parse(merged_file,"fastq")
  variants = [] 
  read_count = 0
  for record in records:
    read_count += 1
    seq  = record.seq
    if len(seq) != 471: continue
    trim_seq  = seq[trim_5:-trim_3]
    if 'N' in trim_seq: continue
    pep = translation(trim_seq)
    variants.append(pep)
  return Counter(variants)

def SD93mut2ID(mut,refseq):
  shift = 248
  haplo = []
  assert(len(mut)==len(refseq))
  for n in range(len(mut)):
    pos = n+shift
    if refseq[n]!=mut[n]:
       haplo.append(refseq[n]+str(pos)+mut[n])
  return '-'.join(haplo)

def Output(InputDict, Rep1Dict, Rep2Dict, outfile, refseq):
  print ("Compiling results into %s" % outfile)
  outfile = open(outfile,'w')
  muts = list(set(list(InputDict.keys())+list(Rep1Dict.keys())+list(Rep2Dict.keys())))
  Rep1_total_count = sum(Rep1Dict.values())
  Rep2_total_count = sum(Rep2Dict.values())
  outfile.write("\t".join(['mut','mutclass','InputCount',
                           'Rep1Count','Rep2Count',
                           'Rep1Freq','Rep2Freq'])+"\n")
  for mut in muts:
    Rep1Freq = float(Rep1Dict[mut])/float(Rep1_total_count)
    Rep2Freq = float(Rep2Dict[mut])/float(Rep2_total_count)
    R1Rep2Enrich = float(Rep2Dict[mut]+1)/float(InputDict[mut]+1)
    mutclass    = 0 if refseq == mut else distance.hamming(refseq,mut)
    ID          = 'WT' if refseq == mut else SD93mut2ID(mut,refseq)
    outfile.write("\t".join(map(str,[ID,mutclass,InputDict[mut],
                                     Rep1Dict[mut],Rep2Dict[mut],
                                     Rep1Freq,Rep2Freq]))+"\n")
  outfile.close()

def main():
  reffile = 'Fasta/SD93_mutlib_ref.fa' 
  outfile = 'result/SD93_MultiMutLib.tsv'
  trim_5  = 25
  trim_3  = 23
  refseq  = next(SeqIO.parse(reffile,"fasta")).seq
  SD93_InputDict   = ProcessMultilib('fastq/SD93_ipt.assembled.fastq', trim_5, trim_3)
  SD93_Rep1Dict    = ProcessMultilib('fastq/SD93_rep1.assembled.fastq', trim_5, trim_3)
  SD93_Rep2Dict    = ProcessMultilib('fastq/SD93_rep2.assembled.fastq', trim_5, trim_3)
  Output(SD93_InputDict,SD93_Rep1Dict,SD93_Rep2Dict,outfile,refseq)

if __name__ == "__main__":
  main()
