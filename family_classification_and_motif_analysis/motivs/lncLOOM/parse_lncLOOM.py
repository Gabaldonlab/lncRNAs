import sys
import glob
import numpy as np

def calculate_real():
	total_families = 0
	signif_fam = 0
	length = []

	with open(sys.argv[1]) as file_list:
			
		for lncLOOM_file in file_list:
			lncLOOM_file = lncLOOM_file.rstrip()
			total_families += 1
			with open(lncLOOM_file, "r+") as input_file:
				f = input_file.readlines()
				#print(f)
				
				### get the number of original seqs
				
				n_seqs = 0 
				index = 0
				for l in f:
					index += 1
					if "Sequences with extended 5' region (5' Graph Calculation):" in l:
						indeq_seq = index
						
				for el in range(9,indeq_seq+1):
					if f[el].startswith(">"):
						n_seqs += 1
				
				#print(n_seqs,lncLOOM_file)

				idx = -1
				motif_list = []
				
				for l in f:
					spp_info = []
					idx += 1
					if ("Motif " in l) and ("Depth:" in l) and not ("Neighborhood" in l):	
						motif_info = []
						motif_info.append(f[idx])			
						for el in range(idx+1,1000000):
							if f[el].startswith("Motif"):
								#print(f[el])
								break
							elif f[el].startswith("E(i)"):
								stats = f[el].rstrip().replace("E(i)-value=","").replace("P(i)-value=","").replace("E(r)-value=","").replace("   ","").split(" ")
								#print(stats)
								stats_int = [float(i) for i in stats]
								motif_info.append(stats_int)
							elif f[el].startswith(">"):
								spp = f[el].replace(">","")
								#print(spp)
								#print(spp.split(" ")[10])
								spp = " ".join(spp.split())
								spp_info.append(spp.split(" ")[0]+"|"+spp.split(" ")[3])
							else:
								continue
						motif_info.append(spp_info)	
						motif_list.append(motif_info)
								

			for el in motif_list:
				if len(el[2]) == n_seqs:
					#print(el)
					if all(float(i) < 0.05 for i in el[1]):
						signif_fam += 1
						print(lncLOOM_file.split("/")[0])
						#print(el[2][0].split("|")[-1])
						length.append(int(el[2][0].split("|")[-1]))
						break
					#print(motif_info)
					#print(f[motif[0]],f[motif[0]+2],f[motif[0]+6],f[motif[1]])
					
					
	#print(signif_fam,"\t",float(signif_fam)/float(total_families)*100)
	#print(length)
	#print(min(length),max(length))
#calculate_real()



def calculate_simulation():
	all_signif_fam = []
	for i in range(1,101):
		loom_files =  glob.glob("simulations/iter%s_*/LincMotif/lncLOOM_Results.txt"%(i))
		total_families = 0
		signif_fam = 0
		length = []

		for lncLOOM_file in loom_files:
			total_families += 1
			with open(lncLOOM_file, "r+") as input_file:
				f = input_file.readlines()
				
				n_seqs = 0 
				index = 0
				for l in f:
					index += 1
					if "Sequences with extended 5' region (5' Graph Calculation):" in l:
						indeq_seq = index
						
				for el in range(9,indeq_seq+1):
					if f[el].startswith(">"):
						n_seqs += 1
						
				idx = -1
				motif_list = []
				
				for l in f:
					spp_info = []
					idx += 1
					if ("Motif " in l) and ("Depth:" in l) and not ("Neighborhood" in l):	
						motif_info = []
						motif_info.append(f[idx])			
						for el in range(idx+1,1000000):
							if f[el].startswith("Motif"):
								#print(f[el])
								break
							elif f[el].startswith("E(i)"):
								stats = f[el].rstrip().replace("E(i)-value=","").replace("P(i)-value=","").replace("E(r)-value=","").replace("   ","").split(" ")
								#print(stats)
								stats_int = [float(i) for i in stats]
								motif_info.append(stats_int)
							elif f[el].startswith(">"):
								spp = f[el].replace(">","")
								#print(spp)
								#print(spp.split(" ")[10])
								spp = " ".join(spp.split())
								spp_info.append(spp.split(" ")[0]+"|"+spp.split(" ")[3])
							else:
								continue
						motif_info.append(spp_info)	
						motif_list.append(motif_info)
						
			for el in motif_list:
				if len(el[2]) == n_seqs:
					#print(el)
					if all(float(s) < 0.05 for s in el[1]):
						signif_fam += 1
						#print(lncLOOM_file,el)
						#print(el[2][0].split("|")[-1])
						#length.append(int(el[2][0].split("|")[-1]))
						break
					#print(motif_info)
					#print(f[motif[0]],f[motif[0]+2],f[motif[0]+6],f[motif[1]])
		print(signif_fam, i)
		all_signif_fam.append(signif_fam)
			
	print(all_signif_fam)
	print(float(np.mean(all_signif_fam))/423*100)
	for el in all_signif_fam:
		if int(el) > 85:
			print(el)
calculate_simulation()
		
