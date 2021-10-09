### This script will generate strand_info.txt file for each species. It has to be run when the pseudomapping step is finished.


import sys
import os


with open("file_list.txt", "r+") as file_list, open("strand_info.txt", "w") as output:
	
	output.write("%s\t%s\t%s\t%s\n"%("sample", "mapping_rate", "strandedness","note"))
	
	for line in file_list:
		line=line.rstrip()

		
		with open(line+"/logs/salmon_quant.log","r+") as log_file:
			strand_info=0
			for line2 in log_file.readlines():
				if "Mapping rate =" in line2:
					mapping_rate = line2.rstrip().split(" ")[-1]
				elif  "Automatically detected most" in line2:
					strand_info+=1
					strand = line2.rstrip().split(" ")[-1]
					

			if strand_info == 0:
				output.write("%s\t%s\t%s\t%s\n"%(line, mapping_rate, "NA", "not used for transcript reconstruction"))
			elif "U" in strand:
				output.write("%s\t%s\t%s\t%s\n"%(line, mapping_rate, strand, "not used for transcript reconstruction"))
			else:
				output.write("%s\t%s\t%s\t%s\n"%(line, mapping_rate, strand, "used for transcript reconstruction"))
			
