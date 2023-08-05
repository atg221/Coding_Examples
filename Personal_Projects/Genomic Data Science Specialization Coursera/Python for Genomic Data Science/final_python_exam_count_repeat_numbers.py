#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:32:17 2023

@author: dg
"""
from collections import Counter
dna=open('/home/dg/Downloads/dna2.fasta')
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

    
#find the number of the repeat in question in each sequence of the fasta file
repeats={}
for i in seqs.keys():
    repeats[i]=[]
    sub_list = []
    for j in range(len(seqs[i]) - 3 ):
        sub_list.append(seqs[i][j:(j+3)])
    print(i)
    print(Counter(sub_list).most_common(5))
    print()
    repeats[i].append(sub_list)  
dna.close()