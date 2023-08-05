#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 14:51:42 2023

@author: dg
"""
dna=open('/home/dg/Downloads/dna2.fasta')
seqs={}
lengths = []
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
#put the lengths of all the sequences in a list
for i in seqs.keys():
    lengths.append(len(seqs[i]))
#sort list in reverse to get largest first
lengths.sort(reverse=True)
print(lengths)

dna.close()