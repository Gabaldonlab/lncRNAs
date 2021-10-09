import sys
import os
import re


with open(sys.argv[1], "r+") as gtf, open(sys.argv[2], "r+") as non_coding_list, open(sys.argv[3], "w") as processed_gtf:


	id_table=[]
	
	non_coding_ids={}
	for line in non_coding_list:
		line=line.rstrip()
		#print line
		xloc=re.findall(r"XLOC_[0-9]{6}", line)
		#print xloc
		non_coding_ids[xloc[0]]=line
	#print non_coding_ids
	
	
	for line in gtf:
		if line.startswith("#"):
			continue
		line=line.rstrip().split("\t")
		if line[2]=="transcript":
			col9=line[8].split(";")
			class_code = col9[-3].split(" ")[-1][1:-1]
			if class_code == "=":
				processed_gtf.write("%s\n"%("\t".join(line)))
			elif class_code == "u":
				processed_gtf.write("%s\n"%("\t".join(line)))
			elif class_code == "x":
				processed_gtf.write("%s\n"%("\t".join(line)))
