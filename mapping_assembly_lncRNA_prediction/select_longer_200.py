import sys
import os
from Bio import SeqIO
import re


with open(sys.argv[1], "r+") as input_fasta, open(sys.argv[2], "w") as output_fasta:
	seqs=[]
	for record in SeqIO.parse(input_fasta, "fasta"):
		if len(record) > 200:
			#print record.id, len(record)
			record.description = record.id #### this outputs only record.id
			seqs.append(record)
			
	SeqIO.write(seqs,output_fasta,"fasta")
