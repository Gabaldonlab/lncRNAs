from Bio import SeqIO
from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA
import sys
import os


retained_sequences=[]
N_retained_seqs=0

discarded_sequences=[]
N_discarded_seqs=0

with open(sys.argv[1], "r") as input_fasta, open(sys.argv[2], "w") as output_fasta: 
	for record in SeqIO.parse(input_fasta, "fasta"):
		if set(record.seq.upper()) <= set(IUPACUnambiguousDNA.letters):
			retained_sequences.append(record)
			N_retained_seqs+=1
		else:
			discarded_sequences.append(record.id)
			N_discarded_seqs+=1 
		
	#print "Discarded sequences"
	
	for el in discarded_sequences:
		print el
	#print N_discarded_seqs
	
	
	#print "Retained sequences"
	#print N_retained_seqs

	SeqIO.write(retained_sequences,output_fasta,"fasta")
