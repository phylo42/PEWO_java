#!/usr/bin/env python

import sys
import re
import collections
print("Conversion from psiblast to fasta")
print("Usage: psiblast2fasta.py aln.psiblast aln.fasta")

f_in=open(sys.argv[1],"r")
f_out=open(sys.argv[2],"w")

#dict to see which header has already been found
headers=collections.OrderedDict()

# read psiblast alignment
lines = f_in.readlines()
for line in lines:
	#skip empty lines
	if (len(line.strip())<1):
		continue;
	#load sequences
	elts=line.strip().split()
	#elts=re.split('\s+',line) //ajout '' before \n, for whatever obsure reason
	if ( elts[0] not in headers ):
		headers[elts[0]]=elts[1]
	else:
		headers[elts[0]]=headers[elts[0]]+elts[1]

#write in output file
for key in headers.keys():
	f_out.write(">"+key+"\n"+headers[key]+"\n")

f_in.close()
f_out.close()

print("DONE!")
