#!/usr/bin/python

def convert_ID_seq(mutlist, muts):
  seq = ''
  for mut_ref in muts:
    pos = mut_ref[1:-1]
    presence = 0
    for mut in mutlist.rsplit('-'):
      if pos in mut:
        seq+=mut[-1]
        presence = 1
        break
    if presence == 0:
      seq+=mut_ref[0]
  return (seq)

def formatting(infile,outfile,cutoff,muts):
  print ("writing: %s" % outfile)
  infile  = open(infile,'r')
  outfile = open(outfile,'w')
  for line in infile.readlines(): 
    if 'mut' in line:
      outfile.write(line.rstrip()+"\t"+"denovo"+"\n")
    else:
      info = line.rstrip().rsplit("\t")
      mutlist  = info[0]
      mutclass = int(info[1])
      iptcount = int(info[2])
      R1count  = int(info[3])
      R2count  = int(info[4])
      R1freq   = float(info[5])
      R2freq   = float(info[6])
      R1enrich = float(info[7])
      R2enrich = float(info[8])
      denovo = []
      if iptcount < cutoff:
        continue 
      for mut in mutlist.rsplit('-'):
        if mut not in muts:
          denovo.append(mut)
      denovo = 'no' if len(denovo) == 0 else '-'.join(denovo)
      if denovo == "no":
        outfile.write(line.rstrip()+"\t"+denovo+"\n")
      if R1enrich > 2 and R2enrich > 2:
        seq = convert_ID_seq(mutlist, muts)
        print (seq)
  infile.close()
  outfile.close()

def main():
  cutoff = 10
  muts    = ['R249K', 'D286G', 'I302V', 'M307V', 'D329N', 'K334S', 'N336H', 'N356D', 'L370S']
  infile  = 'result/Bil69_MultiMutLib.tsv'
  outfile = 'result/Bil69_MultiMutLib_filtered.tsv'
  formatting(infile, outfile, cutoff, muts)

if __name__ == "__main__":
  main()
