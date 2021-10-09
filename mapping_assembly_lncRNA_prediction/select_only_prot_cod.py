import subprocess
from Bio import SeqIO
import sys


record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))
#print(record_dict)


with open(sys.argv[2],"r+") as gene_ids, open(sys.argv[3],"w") as out_fasta:
    gene_ids_list=[]
    
    for l in gene_ids:
        l = l.rstrip()
        gene_ids_list.append(l)
    
    for gene in gene_ids_list:
        for k,v in record_dict.items():
            if gene in k:
                #print(k.split("|")[0:8])
                out_fasta.write(">%s\n%s\n"%("_".join(k.split("|")[0:8]),v.seq))
                
                
                
