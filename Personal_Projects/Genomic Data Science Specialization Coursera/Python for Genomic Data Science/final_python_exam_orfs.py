#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:16:46 2023

@author: dg
"""
dna=open('/home/dg/Downloads/dna2.fasta')
from Bio.Seq import Seq
seqs={}
#make dictionary of IDs and sequences
for line in dna:
    #discard the newline at the end
    line = line.rstrip()
    #distinguish header from sequence
    if line.startswith('>'):
        words=line.split()
        name=words[0][1:]
        seqs[name]=''
    else:
        seqs[name] = seqs[name] + line

#find the orfs of each sequence using the values from the dictionary and print them
for i in seqs.keys():
    print(i)
    seq = Seq(seqs[i])
    table = 1
    min_pro_len = 1
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            for pro in nuc[frame:].translate(table).split("*"):
                if len(pro) >= min_pro_len:
                    print ("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame))
    print()
    
#store the calculated orfs in a new dictionary and the lengths specifically in a seperate one
orfs={}
orflen={}
for i in seqs.keys():
    orfs[i]=[]
    orflen[i]=[]
    seq = Seq(seqs[i])
    table = 1
    min_pro_len = 1
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            for pro in nuc[frame:].translate(table).split("*"):
                if len(pro) >= min_pro_len:
                    orfs[i].append("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame))
                    orflen[i].append(len(pro))

#sort the lengths stored in each dictionary value largest to smallest
for i in orflen.keys():
    orflen[i].sort(reverse=True)
print(orflen)

dna.close()