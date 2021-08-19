# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import sys
from Bio import SeqIO
import subprocess
import random
import numpy as np
import re


def generate_u_fasta(path):
    out_dict={}
    in_dict=SeqIO.to_dict(SeqIO.parse(path,"fasta"))
    for tr_id,seq in in_dict.items():
        if "|u|" in tr_id:
            out_dict[tr_id]=seq
    return(out_dict)


### you can increase the number of simulations by increasing n_iter below (in the paper we used 100 simulations)
n_iter=1+2
def generate_fasta():
    subprocess.call("rm ./*fasta",shell=True)
    family_dict={}
    with open("./intergenic.fam","r+") as families:
        for line in families:
            line=line.rstrip().split("\t")
            #print(line)
            if not line[0] in family_dict:
                family_dict[line[0]] = [line[1]]
            else:
                family_dict[line[0]].append(line[1])
                
    #for k,v in family_dict.items():
    #   print(k,v)
    

    
    calb_u_dict=generate_u_fasta("../blast/calb_lncRNAs.fasta")
    cglab_u_dict=generate_u_fasta("../blast/cglab_lncRNAs.fasta")
    cpar_u_dict=generate_u_fasta("../blast/cpar_lncRNAs.fasta")
    caur_u_dict=generate_u_fasta("../blast/caur_lncRNAs.fasta")
    ctrop_u_dict=generate_u_fasta("../blast/ctrop_lncRNAs.fasta")
	

#    for k,v in ctrop_u_dict.items():
#        print(k,v)
    
    for fam,seqs in family_dict.items():
        with open("./real/%s.fasta"%(fam),"a+") as fam_fasta:
            for seq in seqs:
                if seq in calb_u_dict:
                    fam_fasta.write(calb_u_dict[seq].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                    #print(calb_dict[seq].format("fasta").replace("MSTRG","calb_MSTRG"))
                elif seq in cglab_u_dict:
                    fam_fasta.write(cglab_u_dict[seq].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                elif seq in cpar_u_dict:
                    fam_fasta.write(cpar_u_dict[seq].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                elif seq in caur_u_dict:
                    fam_fasta.write(caur_u_dict[seq].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                else:
                    fam_fasta.write(ctrop_u_dict[seq].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                    

    ### Simulated families
    for iteration in range(1,n_iter):
        print(iteration)
        if not os.path.exists("./simulations/iter_%s"%(iteration)):
			os.makedirs("./simulations/iter_%s"%(iteration))

        for fam,seqs in family_dict.items():
            with open("./simulations/iter_%s/iter%s_random_%s.fasta"%(iteration,iteration,fam),"a+") as random_fam_fasta:
                for seq in seqs:
                    if seq.split("|")[-1] == "calb":
						new_dict = {k: v for k, v in calb_u_dict.items() if k != seq}
						random_fam_fasta.write(new_dict[random.choice(new_dict.keys())].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                    elif seq.split("|")[-1] == "cglab":
						new_dict = {k: v for k, v in cglab_u_dict.items() if k != seq}
						random_fam_fasta.write(cglab_u_dict[random.choice(cglab_u_dict.keys())].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                    elif seq.split("|")[-1] == "cpar":
						new_dict = {k: v for k, v in cpar_u_dict.items() if k != seq}
						random_fam_fasta.write(cpar_u_dict[random.choice(cpar_u_dict.keys())].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                    elif seq.split("|")[-1] == "caur":
						new_dict = {k: v for k, v in caur_u_dict.items() if k != seq}
						random_fam_fasta.write(caur_u_dict[random.choice(caur_u_dict.keys())].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))
                    else:
						new_dict = {k: v for k, v in ctrop_u_dict.items() if k != seq}
						random_fam_fasta.write(ctrop_u_dict[random.choice(ctrop_u_dict.keys())].format("fasta").replace("MSTRG",seq.split("|")[-1]+"_"+"MSTRG"))


generate_fasta()
                    
                    
def run_meme():
	
    meme_real="for fam in real/*fam*.fasta; do echo ${fam} && meme ${fam} -oc $(echo ${fam}| cut -f 1 -d '.')'_res' -nmotifs 1 -mod oops -revcomp -evt 0.05 -dna; done 2>> output_all"
    subprocess.call(meme_real,shell=True)
    
    meme_sim="for fam in simulations/iter*/*fam*.fasta; do echo ${fam} && meme ${fam} -oc $(echo ${fam}| cut -f 1 -d '.')'_res' -nmotifs 1 -mod oops -revcomp -evt 0.05 -dna; done 2>> output_all"
    subprocess.call(meme_sim,shell=True)
    
    
    
run_meme()


def calculate_e_value():
    
    N_fams=423
    
    dir_list=[]
    for d in os.walk("./real/"):
        if "fam" in d[0]:
            dir_list.append(d[0])


    N_fam_with_motifs_in_real=0
    motiv_list=[]
    for d in dir_list:
        if not "random" in d:
#            print(d[0])
            with open("%s/meme.txt"%(d),"r+") as meme_res:
                lines_list = meme_res.read().splitlines()
                #print(lines_list)
                for line in lines_list:
                    #print(line)
                    #line=line.rstrip()
                    if "E-value =" in line:
                        if float(line.split("=")[-1])<0.05:
                            #print(d,line)
                            cons = [x for x in lines_list if "Multilevel" in x]
                            cons_seq=cons[0].split(" ")[-1]
                            cons_seq=cons_seq.replace("T","U")
                            #print(lines_list[index+2])
                            motiv_list.append([d,"\t",cons_seq])
                            N_fam_with_motifs_in_real+=1
    prcnt_real=float(N_fam_with_motifs_in_real)/float(N_fams)*100
    
    with open("motivs_in_families.txt","w") as output:
        for l in motiv_list:
			print(l[0].split("/")[-1].replace("_res",""))
			output.write("%s\n"%(l[0].split("/")[-1].replace("_res","")))
    

    N_fam_with_motifs_in_random=[]

    dir_list_sim=[]
    for d in os.walk("./simulations/"):
        if "fam" in d[0]:
            dir_list_sim.append(d[0])

    for i in range(1,n_iter):
        #print(i)
        N_fam_with_motifs_in_rand=0
        iter_list = [k for k in dir_list_sim if "iter%s_"%(i) in k]
        for d in iter_list:
            #print(d)
            with open("%s/meme.txt"%(d),"r+") as meme_res:
                for line in meme_res:
                    line=line.rstrip()
                    if "E-value =" in line:
                        if float(line.split("=")[-1])<0.05:
                            print(d,line)
                            N_fam_with_motifs_in_rand+=1
        prcnt_random=float(N_fam_with_motifs_in_rand)/float(N_fams)*100
        N_fam_with_motifs_in_random.append(prcnt_random)
    print(N_fam_with_motifs_in_random)

    n=0
    for val in N_fam_with_motifs_in_random:
		if val>=prcnt_real:
			n+=1
			
    #print(np.mean(N_fam_with_motifs_in_random))
    with open("output_final_results_test.txt","w") as final_res:
		final_res.write("%s simulation sets are more/equal than %s%%, pvalue = %s"%(n,round(prcnt_real,3),float(n)/n_iter))
		final_res.write("real data = %s percent, simulated = %s (sd=%s) percent"%(round(prcnt_real,3),round(np.mean(N_fam_with_motifs_in_random),3),round(np.std(N_fam_with_motifs_in_random),3)))

		print("%s simulation sets are more/equal than %s%%, pvalue = %s"%(n,round(prcnt_real,3),float(n)/n_iter))
		print("real data = %s percent, simulated = %s (sd=%s) percent"%(round(prcnt_real,3),round(np.mean(N_fam_with_motifs_in_random),3),round(np.std(N_fam_with_motifs_in_random),3)))

calculate_e_value()
