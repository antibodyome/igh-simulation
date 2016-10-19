import sys
from Bio import Seq,SeqIO

## Set parameters

# Outfile name
oname=str(sys.argv[1])

## Shouldn't have to change anything below

## Read in V, D, and J

vname="IGHV.fas"
dname="IGHD.fas"
jname="IGHJ.fas"

vseqs={}
for record in SeqIO.parse(vname,"fasta"):
  vname=record.id
  vallele=vname.split("*")[1]
  if vallele=="01":
    vseqs[vname]=str(record.seq).upper()
vgenes=list(vseqs.keys())
vgenes.sort()

dseqs={}
for record in SeqIO.parse(dname,"fasta"):
  dname=record.id
  dallele=dname.split("*")[1]
  if dallele=="01":
    dseqs[record.id]=str(record.seq).upper()
dgenes=list(dseqs.keys())
dgenes.sort()

jseqs={}
for record in SeqIO.parse(jname,"fasta"):
  jname=record.id
  jallele=jname.split("*")[1]
  if jallele=="01":
    jseqs[record.id]=str(record.seq).upper()
jgenes=list(jseqs.keys())
jgenes.sort()

ofile=open(oname,'w')
i=1
for vg in vgenes:
  v=vseqs[vg]
  for dg in dgenes:
    d=dseqs[dg]
    for jg in jgenes:
      j=jseqs[jg]
      ofile.write(">S"+str(i)+"|"+vg+"|"+dg+"|"+jg+"\n"+v+d+j+"\n")
      i+=1
ofile.close()
