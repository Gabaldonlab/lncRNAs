#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 12:00:31 2020

@author: hhovhannisyan
"""


import subprocess
import sys
import os
import itertools
from Bio import SeqIO


species_u=["calb_u","ctrop_u","cpar_u","caur_u","cglab_u"]

### Run RNAfold
for spp in species_u:
    rnafold_u="RNAfold --noPS -i ./%s_lncRNAs_renamed.fasta -j10 > %s_lncRNAs_fold.fastb"%(spp,spp)
    subprocess.call(rnafold_u,shell=True)

### remove extra character from files
remove_brackets='for fold in *fastb; do cut -f1 -d " " ${fold} > final_${fold};done'
subprocess.call(remove_brackets,shell=True)


### parse RNAfold files to fix Beagle bug
def fix_RNAfold(species):
    for spp in species:
        for spp_file in [spp]:
            spp_dict={}
            with open("final_"+spp_file+"_lncRNAs_fold.fastb","r+") as fold_file, open("final_"+spp_file+"_lncRNAs_fold_for_beagle.fastb","w") as output:
                #spp_file+"_lncRNAs_fold.fastb
                lines=fold_file.readlines()
                #print(lines)
                index=0
                while index < len(lines):
                    spp_dict[lines[index]]=[lines[index+1],lines[index+2]]
                    #print(index,index+1,index+2)
                    index+=3
                    
                for seq_id,seq_struct in spp_dict.items():
                    if seq_struct[1].startswith("("):
                        output.write("%s%s%s"%(seq_id,"N"+seq_struct[0],"."+seq_struct[1]))
                        #print(seq_id)
                        #print("N"+seq_struct[0])
                        #print("."+seq_struct[1])
                    else:
                        output.write("%s%s%s"%(seq_id,seq_struct[0],seq_struct[1]))

fix_RNAfold(species_u)
                    
                    
###### run beagle in pairwise manner. This will create batch_beagle.txt file with all pairwise commands. We recommend to use HPC to run it.
all_beagle_runs=[]
def run_pairwise_beagle(spp_list):

    all_pairwise_comparisons_without_repeats_and_without_doubles=list(itertools.combinations(spp_list, 2))
    print(all_pairwise_comparisons_without_repeats_and_without_doubles)
    batch_beagle_commands=[]
    
    
    for comparison in all_pairwise_comparisons_without_repeats_and_without_doubles:
        spp1=comparison[0]
        spp2=comparison[1]
        if not os.path.exists('./%s_%s/'%(spp1,spp2)):
			os.makedirs('./%s_%s/'%(spp1,spp2))

        with open("final_%s_lncRNAs_fold_for_beagle.fastb"%(spp1),"r+") as spp1_file:
        #with open("test.fastb") as spp1_file:
			spp1_file = spp1_file.readlines()
			n_rows = len(spp1_file)

			
			start=0
			
			for i in range(1,(n_rows/3)+1):
				with open("%s_%s/%s_%s.fastb"%(spp1,spp2,spp1,i),"w") as out_f:
					to_write = spp1_file[start:start+3]
					for el in to_write:
						out_f.write(el)
					start+=3
					batch_beagle_commands.append("bash %s_%s/%s_%s_to_%s.sh"%(spp1,spp2,spp1,i,spp2))
					with open("%s_%s/%s_%s_to_%s.sh"%(spp1,spp2,spp1,i,spp2),"w") as beagle_command:
						beagle_command.write("java -jar -Xloggc:%s_%s/%s_%s_to_%s.log Beagle_v0.2.jar -input1 %s_%s/%s_%s.fastb -input2 final_%s_lncRNAs_fold_for_beagle.fastb -c 2 -l true -pValue true -outfile %s_%s/%s_%s_to_%s_beagle_outfile.txt\n"%(spp1,spp2,spp1,i,spp2,spp1,spp2,spp1,i,spp2,spp1,spp2,spp1,i,spp2))
					
	if "_u" in spp1:
		with open("batch_beagle_u_commands.txt","w") as out_beagle:
			for el in batch_beagle_commands:
				out_beagle.write("%s\n"%(el))

run_pairwise_beagle(species_u)


###### grep ids from each file and make a concatenated file

def grep_ids(spp_list):
	
	if os.path.exists("all_beagle_results_ids.txt"):
		os.remove("all_beagle_results_ids.txt")
		
	if os.path.exists("all_beagle_results_ids_u.txt"):
		os.remove("all_beagle_results_ids_u.txt")
		
	all_pairwise_comparisons_without_repeats_and_without_doubles=list(itertools.combinations(spp_list, 2))
	print(all_pairwise_comparisons_without_repeats_and_without_doubles)
	
	for comparison in all_pairwise_comparisons_without_repeats_and_without_doubles:
		spp1=comparison[0]
		spp2=comparison[1]
		if "_u" in spp1:
			grep= 'grep ">MSTRG" %s_%s/%s_*_to_%s_beagle_outfile.txt | cut -f 2,3,4,5,6,7,8 -d ":" >> all_beagle_results_ids_u.txt'%(spp1,spp2,spp1,spp2)
			print(grep)
		#run
		subprocess.call(grep,shell=True)
		
	
	if "calb_u" in spp_list:
		file_for_cluster_script = "python parse_beagle.py all_beagle_results_ids_u.txt intergenic . "
	subprocess.call(file_for_cluster_script,shell=True)
	
	

grep_ids(species_u)

        
###### Generate a file for clustering script

def make_family_command(alignments,spp_list,directory):
    N_seqs_list=[]
    for spp in spp_list:
        fastb_file=spp+"_lncRNAs_fold.fastb"
        
        N_seqs=0
        with open(fastb_file, "r+") as f_file:
			for line in f_file:
				N_seqs+=1
        N_seqs_list.append(str(int(N_seqs)/3))
    command="python ./classifyFamiliesv5_VennGH_mod_for_candidas_5spp.py "+ str(" ".join(N_seqs_list))+ " %s "%(alignments) + "%s/%s.fam "%(directory,directory)+"%s/%s.txt "%(directory,directory)+ "%s/%s_FAM_VENN.R "%(directory,directory)+"%s/%s_GENES_VENN.R "%(directory,directory)+"> %s/output_%s_families "%(directory,directory)
    
    print(command)
    subprocess.call(command,shell=True)

make_family_command("results_for_clustering_script_intergenic.txt",species_u, "intergenic")




### run r scripts to make venn diagramms
run_r_scripts = "Rscript intergenic/intergenic_FAM_VENN.R && \
    Rscript intergenic/intergenic_GENES_VENN.R"
subprocess.call(run_r_scripts,shell=True)   

