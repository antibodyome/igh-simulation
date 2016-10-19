# This checks for:
# CDR3 defined by a regex
# no stops,
# the TVSS motif at the end
# the conserved J motif

import sys
import csv
import re
import numpy
from collections import Counter
from Bio import Seq,SeqIO
from progress.bar import Bar

## Set parameters

# Number of reads
numreads=int(sys.argv[1])

# Outfile name stub
onamestub=str(sys.argv[2])

# Seed
seed=int(sys.argv[3])

## Shouldn't have to change anything below

# Set regex for CDR3; make sure it is present in all reads
cdr3p=re.compile("(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]")

# J region motifs
jregionre=re.compile("[FW]G[A-Z]G")
#jendre=re.compile("T[LMT]VTVSS")

# Numerical stuff
numpy.random.seed(seed)
rearrange_attempts=10
insert_attempts=10

# Base frequencies for insertion
# aken from Jackson et al 2004
bases=['A','C','G','T']
inscount=[530,803,998,551]
insfreq=[float(b)/sum(inscount) for b in inscount]

## Read in V, D, and J

vname="IGHV.fas"
dname="IGHD.fas"
jname="IGHJ.fas"

vseqs={}
for record in SeqIO.parse(vname,"fasta"):
  vseqs[record.id]=str(record.seq).upper()
vgenes=list(vseqs.keys())
vgenes.sort()
vgenecount=Counter([s.split("*")[0] for s in vgenes])
vgenep=[]
for i in range(len(vgenes)):
  vg=vgenes[i]
  vgc=vgenecount[vg.split("*")[0]]
  vgenep.append(1.0/vgc)
vgenep=[x/sum(vgenep) for x in vgenep]


dseqs={}
for record in SeqIO.parse(dname,"fasta"):
  dseqs[record.id]=str(record.seq).upper()
dgenes=list(dseqs.keys())
dgenes.sort()
dgenecount=Counter([s.split("*")[0] for s in dgenes])
dgenep=[]
for i in range(len(dgenes)):
  dg=dgenes[i]
  dgc=dgenecount[dg.split("*")[0]]
  dgenep.append(1.0/dgc)
dgenep=[x/sum(dgenep) for x in dgenep]


jseqs={}
for record in SeqIO.parse(jname,"fasta"):
  jseqs[record.id]=str(record.seq).upper()
jgenes=list(jseqs.keys())
jgenes.sort()
jgenecount=Counter([s.split("*")[0] for s in jgenes])
jgenep=[]
for i in range(len(jgenes)):
  jg=jgenes[i]
  jgc=jgenecount[jg.split("*")[0]]
  jgenep.append(1.0/jgc)
jgenep=[x/sum(jgenep) for x in jgenep]

# Read in VDJ file from Jackson et al.

vdjfile=open("vdj.csv",'r')

v3del=[]
d3del=[]
d5del=[]
j5del=[]
vdins=[]
djins=[]

vdjreader=csv.reader(vdjfile,delimiter=",")
vdjheader=next(vdjreader)
for row in vdjreader:
  v3del.append(int(row[11]))
  d3del.append(int(row[13]))
  d5del.append(int(row[14]))
  j5del.append(int(row[12]))
  vdins.append(int(row[9]))
  djins.append(int(row[10]))

reads=[]
vdj=[]
labels=[]
bar = Bar('Generating', max=numreads)
for read in range(numreads):
  rattempt=0
  while True:
    ## Rearrangements
    vg=numpy.random.choice(vgenes,None,False,vgenep)
    v=vseqs[vg]
    dg=numpy.random.choice(dgenes,None,False,dgenep)
    d=dseqs[dg]
    jg=numpy.random.choice(jgenes,None,False,jgenep)
    j=jseqs[jg]

    ## Indels
    iattempt=0
    while True:
      # Deletions
      v3d=numpy.random.choice(v3del)
      d5d=numpy.random.choice(d5del)
      d3d=numpy.random.choice(d3del)
      j5d=numpy.random.choice(j5del)

      # Insertions
      vdi=numpy.random.choice(vdins)
      dji=numpy.random.choice(djins)

      ## Make sequence
      # First delete 3' of V
      vs=v[:len(v)-v3d]
      # ..both ends of D
      ds=d[d5d:len(d)-d3d]
      # ...and 5' of J
      js=j[j5d:]

      # Insertions
      vdis=''.join(numpy.random.choice(bases,vdi,True,insfreq))
      djis=''.join(numpy.random.choice(bases,dji,True,insfreq))

      # Now join
      igh=vs+vdis+ds+djis+js
      label=[len(vs),len(vdis),len(ds),len(djis),len(js)]
      # Check if rearranged sequence is coding
      ighp=Seq.Seq(igh).translate()
      numstops=ighp.count('*')
      # Check if sequence has a CDR3
      m=cdr3p.search(igh)
      # Check that the sequence ends with TVSS
      ighptail=str(ighp)[-4:]
      # Check for conserved J motif
      m2=jregionre.search(str(ighp))
      if numstops==0 and m!=None and m2!=None and ighptail=="TVSS":
        break
      iattempt+=1
      if iattempt==insert_attempts:
        break
    if numstops==0 and m!=None and m2!=None and ighptail=="TVSS":
      reads.append(igh)
      vdj.append([vg,dg,jg])
      labels.append(label)
      bar.next()
      break
    else:
      rattempt+=1
    if rattempt==rearrange_attempts:
      print("Failed rearrangement")
      break

bar.finish()
totalreads=len(reads)
seqnames=["S"+str(i+1) for i in range(totalreads)]
ofile=open(onamestub+".fas",'w')
for i in range(totalreads):
  germline=vdj[i]
  read=reads[i]
  ofile.write(">"+seqnames[i]+"|"+germline[0]+"|"+germline[1]+"|"+germline[2]+"\n"+read+"\n")
ofile.close()

ofile=open(onamestub+".txt",'w')
ofile.write("seqname\tV\tD\tJ\tvlen\tvdilen\tdlen\tdjilen\tjlen\n")
for i in range(totalreads):
  germline=vdj[i]
  regions=[str(l) for l in labels[i]]
  ofile.write(seqnames[i]+"\t"+germline[0]+"\t"+
    germline[1]+"\t"+germline[2]+"\t"+
    regions[0]+"\t"+regions[1]+"\t"+regions[2]+"\t"+regions[3]+"\t"+regions[4]+"\n")
ofile.close()
