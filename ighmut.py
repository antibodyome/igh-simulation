import sys
import csv
import numpy
from Bio import Seq,SeqIO
from progress.bar import Bar

## Set parameters

# Number of mutations
nummut=int(sys.argv[1])

# I/O
iname=str(sys.argv[2])
oname=str(sys.argv[3])

# Seed
seed=int(sys.argv[4])

## Shouldn't have to change anything below
# Set seed
numpy.random.seed(seed)

# Attempts
mutate_attempts=100

## Read in S5F model

mr=csv.reader(open("Mutability.csv",'r'), delimiter=' ')
mrheader=next(mr)
global mutdict
mutdict={}
for row in mr:
  mutdict[row[0]]=float(row[1])

# This refers to the columns of the substitution matrix
bases=['A','C','G','T']
sr=csv.reader(open("Substitution.csv",'r'), delimiter=' ')
srheader=next(sr)
global substdict
substdict={}
for row in sr:
  substdict[row[0]]=[float(row[1]),float(row[2]),float(row[3]),float(row[4])]

## Read in unmutated sequences

seqs=[]
for record in SeqIO.parse(iname,"fasta"):
  seqs.append(record)
numseq=len(seqs)

def get_fivemers(s):
  fivemers=[]
  for i in range(len(s)-4):
    fivemers.append(s[i:i+5])
  return fivemers

def get_mutability(fivemers):
  mut=[mutdict[f] for f in fivemers]
  return [0,0]+mut+[0,0]

## Main loop
# Initialise array of mutated sequences and IDs
mutseqs=[]
mutseqids=[]
bar = Bar('Mutating', max=numseq)
for s in range(numseq):
  strseq=str(seqs[s].seq)
  strseqchars=[c for c in strseq]
  lseq=len(strseq)
  mutseq=strseq
  m=0
  while m<nummut:
    # Initialise mutated sequence
    # Get fivemers and mutability
    fivemers=get_fivemers(mutseq)
    mut=get_mutability(fivemers)
    # Normalise mutability
    normmut=[m/sum(mut) for m in mut]
    # Try mutating the current sequence
    mtry=0
    while True:
      # Identify position to be mutated
      pos=numpy.random.choice(range(lseq),None,replace=True,p=normmut)
      this5mer=fivemers[pos-2]
      substprob=substdict[this5mer]
      nuc=str(numpy.random.choice(bases,None,replace=True,p=substprob))
      # Check new sequence is coding
      newmutseq=''.join([mutseq[:pos],nuc,mutseq[(pos+1):]])
      newmutseqp=Seq.Seq(newmutseq).translate()
      numstops=newmutseqp.count('*')
      if numstops==0:
        mutseq=newmutseq
        break
      mtry+=1
      if mtry==mutate_attempts:
        print("Reached mutate attempts")
        break
    mutseqchars=[c for c in mutseq]
    muts=[mutseqchars[p]!=strseqchars[p] for p in range(len(strseq))]
    m=sum(muts)
  bar.next()
  if mtry<mutate_attempts:
    mutseqs.append(mutseq)
    mutseqids.append(seqs[s].id)

bar.finish()

mutseqnames=[i+"|"+str(nummut) for i in mutseqids]
ofile=open(oname,'w')
for i in range(len(mutseqs)):
  ofile.write(">"+mutseqnames[i]+"\n"+mutseqs[i]+"\n")
ofile.close()
