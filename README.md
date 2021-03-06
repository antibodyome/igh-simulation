# IGH simulation

This repository contains Python scripts to simulate IGH rearrangements, in order to test V(D)J assignment utilities. These were used to generate simulated data for [our paper](http://rstb.royalsocietypublishing.org/content/370/1676/20140240.long) on V(D)J assignment using phylogenetic placement. The code is very simple (and in places, horribly inefficient), and at some point, this code may be overhauled.

It uses the following files:
- ```IGHV.fas```, ```IGHD.fas```, ```IGHJ.fas```: these are functional V, D and J human genes obtained from IMGT (IMGT/V-QUEST reference directory release 201443-5, 24 October 2014). Ambiguous nucleotides have been resolved to the base of the most closely related allele. Some alleles contain Ns.
  - IGHV1-2*03 has a N replaced with G
  - IGHV1-45*01 has a N replaced with T
  - IGHV4-30-4*05 N replaced with A
  - IGHV4-4*06 has NNN replaced with AGG
  - IGHV4-59*09 has NANNN replaced with GAAGG
- ```vdj.csv```: this contains information on inferred V(D)J rearrangements from a curated set of immunoglobulin sequences analysed by [Jackson *et al.* (2004)](http://www.biomedcentral.com/1471-2172/5/19). This was taken from an [Excel sheet](http://www.biomedcentral.com/content/supplementary/1471-2172-5-19-s1.xls) in the supplementary information. I renamed the columns, and saved as a csv file.
- ```Mutability.csv``` and ```Substitution.csv```: these are files on mutability and substitution probabilities from the S5F model of [Yaari et al. 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3828525/), version 0.1 (July 26, 2013), downloaded from [the Kleinstein lab webpage](http://clip.med.yale.edu/shm/download.php).

Dependencies:
- Python
- Python libraries (e.g. installable via `pip`)
  - Biopython
  - Numpy
  - progress

## ```ighshuffle.py```

This Python script generates all possible V(D)J rearrangements of *01 alleles, without insertions, deletions, or mutations.

Usage

```python ighshuffle.py <output filename>```

## ```ighsim.py```

This Python script generates V(D)J rearrangements using simple random sampling (SRS) of the V, D and J regions. Deletion lengths in the V, D and J regions and the lengths of the N1 and N2 non-templated regions are generated by resampling the distributions from Jackson *et al.* (given in the file `vdj.csv`). Base frequencies for N1 and N2 are from concatenated N1 and N2 regions. The resulting sequences are coding (no stops), have the conserved ```[FW]G.G``` motif, end with the J region motif ```TVSS```, and are restricted to have a CDR3 region as defined by the following (Pythonic) regular expression.

    ```(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]```

Usage:

```python ighsim.py <number of reads> <output filename stub> <random number seed>```

This will generate a FASTA file (```<output filename stub>.fas```), as well as a table of germlines (```<output filename stub>.txt```).

## ```ighmut.py```

This script takes a FASTA file, e.g. one produced by ```ighshuffle.py``` or ```ighsim.py```, and adds mutations according to the S5F model (or at least, my interpretation of it). This script ensures that despite mutations, the above conditions (no stops, conserved motifs, and the regular expression) are preserved; there is a hard-coded limit of 100 mutation attempts per sequence.

Usage:

```python ighmut.py <number of mutations> <input fasta> <output fasta> <random number seed>```

## Examples

Some simulated data have been included:

```
python ighshuffle.py simple_rearrangements.fas
python ighsim.py 1000 simple_plus_indels 1
python ighmut.py 40 simple_plus_indels.fas simple_plus_indels_40.fas 1
```

## To do

- Improve efficiency, especially for mutations

## Citation

If for some bizarre reason you use this code in your research, please cite:

Frost SDW, Murrell B, Hossain ASMM, Silverman GJ, Kosakovsky Pond SL. (2015) Assigning and visualizing germline genes in antibody repertoires. *Philos Trans R Soc Lond B Biol Sci.* 370:20140240. [doi:10.1098/rstb.2014.0240](http://dx.doi.org/10.1098/rstb.2014.0240).
