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
      denovo = []
      if max([iptcount, R1count, R2count]) < cutoff:
        continue 
      for mut in mutlist.rsplit('-'):
        if mut not in muts:
          denovo.append(mut)
      denovo = 'no' if len(denovo) == 0 else '-'.join(denovo)
      outfile.write(line.rstrip()+"\t"+denovo+"\n")
      if R1freq > 0.001 and R2freq > 0.001:
        seq = convert_ID_seq(mutlist, muts)
        print (seq)
  infile.close()
  outfile.close()

def main():
  cutoff = 10
  muts    = ['E248G', 'R249K', 'I265T', 'Y336H', 'R338L', 'N339D', 'K344E', 'S346G', 'E369K', 'G381E', 'N387K']
  infile  = 'result/SD93_MultiMutLib.tsv'
  outfile = 'result/SD93_MultiMutLib_filtered.tsv'
  formatting(infile, outfile, cutoff, muts)

if __name__ == "__main__":
  main()
