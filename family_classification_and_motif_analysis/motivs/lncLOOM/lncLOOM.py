import os
import subprocess
import sys
import glob



def run_lncLOOM():
	for fasta_file in glob.glob("../real/*fasta"%(path)):
		#print(os.path.basename(fasta_file))
		fasta = os.path.basename(fasta_file)
		lncLOOM_real = "python LncLOOM.py --fasta ./real/%s --startw 30 -r 100 -o %s -m 8 1>>output 2>>error"%(fasta,os.path.splitext(fasta)[0])
		print(lncLOOM_real)
		subprocess.call(lncLOOM_real, shell=True)
run_lncLOOM()


def run_lncLOOM_sim():
	### you can increase the numer of simualtions (there were 100 in the paper)
	for i in range(1,2): 
		for fasta_file in glob.glob("../simulations/iter%s_*fasta"%(path,i)):
			#print(os.path.basename(fasta_file))
			fasta = os.path.basename(fasta_file)
			lncLOOM_sim = "python LncLOOM.py --fasta ../simulations/%s --startw 30 -r 100 -o ./simulations/iter_%s/%s -m 2 1>>output_sim_mn 2>>error_sim_mn"%(fasta,i,os.path.splitext(fasta)[0])
			print(lncLOOM_sim)
			subprocess.call(lncLOOM_sim, shell=True)
run_lncLOOM_sim()
