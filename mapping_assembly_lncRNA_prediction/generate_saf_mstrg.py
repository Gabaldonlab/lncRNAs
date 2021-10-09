import sys
import os
import re
import numpy


with open(sys.argv[1], "r+") as gtf, open(sys.argv[2], "r+") as non_coding_list, open(sys.argv[4], "w") as saf, open(sys.argv[5], "w") as bed, open(sys.argv[6], "r") as original_gff:

	id_table=[]
	
	non_coding_ids={}
	for line in non_coding_list:
		line_split=line.rstrip().split("_")
		mstrg=line_split[0]
		non_coding_ids[mstrg]=line
	#print non_coding_ids
	
	
	for line in gtf:
		if line.startswith("#"):
			continue
		line=line.rstrip().split("\t")
		if line[2]=="transcript":
			col9=line[8].split(";")
			class_code = col9[-3].split(" ")[-1][1:-1]
			#print class_code
			length=abs(int(line[3])-int(line[4]))+1
			#print length
			if class_code == "=":
				continue


			elif class_code == "u":
				mstrg=col9[0].split(" ")[1][1:-1]
				if mstrg in non_coding_ids:
					mstrg=mstrg+"|"+str(length)+"|"+class_code+"|"+"NA"+"|"+line[0]+"|"+sys.argv[3]
					if not mstrg in id_table:
						gtf_line="%s;%s;%s;%s;%s"%(mstrg, line[0],line[3],line[4],line[6])
						id_table.append(gtf_line)
				

			elif class_code == "x":
				mstrg=col9[0].split(" ")[1][1:-1]
				if mstrg in non_coding_ids:
					gene_name = col9[4].split(" ")[2][1:-1]
					mstrg=mstrg+"|"+str(length)+"|"+class_code+"|"+gene_name+"|"+line[0]+"|"+sys.argv[3]
					if not mstrg in id_table:
						gtf_line="%s;%s;%s;%s;%s"%(mstrg, line[0],line[3],line[4],line[6])
						id_table.append(gtf_line)
						
						
						
				
	for line in original_gff:
		if line.startswith("#"):
			continue
		line=line.rstrip().split("\t")
		if line[2]=="gene":
			ID=line[8].split(";")[0].replace("ID=","")
			#print ID
			orig_gtf_line="%s;%s;%s;%s;%s"%(ID+"|"+str(abs(int(line[3])-int(line[4]))+1), line[0],line[3],line[4],line[6])
			id_table.append(orig_gtf_line)
	
					
	#saf.write("%s\t%s\t%s\t%s\t%s\n"%("GeneID","Chr","Start","End","Strand"))
	for k in id_table:
		k=k.split(";")
		k="\t".join(k)
		saf.write("%s\n"%(k))
		
	#bed.write("%s\t%s\t%s\t%s\t%s\n"%("Chr","Start","End","TransID","Strand"))
	for line in id_table:
		line=line.split(";")
		line_rearranged="%s;%s;%s;%s;%s"%(line[1],line[2],line[3],line[0],line[4])
		line_rearranged_split=line_rearranged.split(";")
		line_rearranged_split="\t".join(line_rearranged_split)
		bed.write("%s\n"%(line_rearranged_split))
		

