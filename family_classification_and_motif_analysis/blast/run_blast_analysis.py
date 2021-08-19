
import subprocess
import sys
import os
import itertools
from Bio import SeqIO


species_u=["calb_u","ctrop_u","cpar_u","caur_u","cglab_u"]


##### make separate files with u lncRNAs
def select_u(spp_list):
    for spp in spp_list:
        with open(spp+"_lncRNAs.fasta","w") as u_lncrnas:
			
			if "ctrop" in spp:
				discard_list = []
				with open("../../ploting_and_network_analysis/ctrop/ctrop_to_discard.txt","r+") as discard_file:
					for l in discard_file:
						l=l.rstrip()
						discard_list.append(l)
			elif "caur" in spp:
				discard_list = []
				with open("../../ploting_and_network_analysis/caur/caur_to_discard.txt","r+") as discard_file:
					for l in discard_file:
						l=l.rstrip()
						discard_list.append(l)
			else:
				print(spp)
				
			if "ctrop" in spp:
				print(discard_list)		
			for seq_record in SeqIO.parse("./%s_lncRNAs.fasta"%(spp.split("_")[0]), "fasta"):
				name = str(seq_record.id)
				seq = str(seq_record.seq)
				if "|u|" in name:
					if "ctrop" in spp:
						if not name in discard_list:
							u_lncrnas.write("%s%s\n"%(">",name))
							u_lncrnas.write("%s\n"%(seq))
					elif "caur" in spp:
						if not name in discard_list:
							u_lncrnas.write("%s%s\n"%(">",name))
							u_lncrnas.write("%s\n"%(seq))
					else:
						u_lncrnas.write("%s%s\n"%(">",name))
						u_lncrnas.write("%s\n"%(seq))
						
						
select_u(species_u)
#####


##### perform makeblastdb

def makeblastbd(spp_list):
    for spp in spp_list:
        rename_ids="sed '/^>/ s/ .*//' %s_lncRNAs.fasta > %s_lncRNAs_renamed.fasta"%(spp,spp)
        print(rename_ids)
        subprocess.call(rename_ids, shell=True)
        
        
        make_blast="makeblastdb -in %s_lncRNAs_renamed.fasta -dbtype nucl -out blastdbs/%s 1>blastdbs/output_%s 2>blastdbs/error_%s "%(spp,spp,spp,spp)
        print(make_blast)
        subprocess.call(make_blast, shell=True)        
        
        

makeblastbd(species_u)
######



###### run blast in pairwise manner

def run_pairwise_blastn(spp_list):

    all_pairwise_comparisons_without_repeats_and_without_doubles=list(itertools.combinations(spp_list, 2))
    print(all_pairwise_comparisons_without_repeats_and_without_doubles)

    for comparison in all_pairwise_comparisons_without_repeats_and_without_doubles:
        spp1=comparison[0]
        spp2=comparison[1]
        ### if there are more than 1 good alignemnt for a query, -max_target_seqs 1 select the one with lowest p-value. 
        ### is object-query pair has more than one alignment, -max_hsps 1 chooses the one with lowest p-value
        run_blast_spp2_to_spp1="blastn -db blastdbs/%s -num_threads 10 -query %s_lncRNAs_renamed.fasta -max_target_seqs 5 -max_hsps 5 -out %s_to_%s.tsv -outfmt 6 -evalue 1e-3"%(spp1,spp2,spp2,spp1)
        subprocess.call(run_blast_spp2_to_spp1,shell=True)
        print(run_blast_spp2_to_spp1)
        
        run_blast_spp1_to_spp2="blastn -db blastdbs/%s -num_threads 10 -query %s_lncRNAs_renamed.fasta -max_target_seqs 5 -max_hsps 5 -out %s_to_%s.tsv -outfmt 6 -evalue 1e-3"%(spp2,spp1,spp1,spp2)
        subprocess.call(run_blast_spp1_to_spp2,shell=True)
        print(run_blast_spp1_to_spp2)
        
    return(all_pairwise_comparisons_without_repeats_and_without_doubles)
        

run_pairwise_blastn(species_u)
########


######## find reciprocal hits


def find_reciprocal_hits(spp_list,outfile):
    
    all_pairwise_comparisons_without_repeats_and_without_doubles=list(itertools.combinations(spp_list, 2))
    with open(outfile,"w") as output_synteny:
        for comparison in all_pairwise_comparisons_without_repeats_and_without_doubles:
        
            spp1_file = comparison[0]+"_to_"+comparison[1]+".tsv"
            spp2_file = comparison[1]+"_to_"+comparison[0]+".tsv"
        
            spp1_hits_list=[]
            spp2_hits_list=[]
            
            #spp1_best_hits_list=[]
            #spp2_best_hits_list=[]
        
            with open(spp1_file,"r+") as spp1_hits, open(spp2_file, "r+") as spp2_hits:
                for line1 in spp1_hits:
                    line1=line1.rstrip().split("\t")
                    hits1=[line1[0],line1[1]]
                    spp1_hits_list.append(hits1)
                
           
                    
                    
                
                for line2 in spp2_hits:
                    line2=line2.rstrip().split("\t")
                    hits2=[line2[0],line2[1]]
                    spp2_hits_list.append(hits2)
                
        
            for el1 in spp1_hits_list:
                for el2 in spp2_hits_list:
                    if set(el1)==set(el2):
                        #print(comparison[0],comparison[1],el1[0],el1[1])
                        output_synteny.write("%s\t%s\t%s\t%s\n"%(comparison[0],comparison[1],el1[0],el1[1]))
                    
            for el2 in spp2_hits_list:
                for el1 in spp1_hits_list:
                    if set(el2)==set(el1):
                        output_synteny.write("%s\t%s\t%s\t%s\n"%(comparison[1],comparison[0],el2[0],el2[1]))
                        #print(comparison[1],comparison[0],el2[0],el2[1])


find_reciprocal_hits(species_u, "intergenic/synteny_info_output_intergenic.txt")

######


###### Generate a file for clustering
commands=[]
def make_family_command(spp_list,directory):
    N_seqs_list=[]
    for spp in spp_list:
        fasta_file=spp+"_lncRNAs.fasta"
        
        N_seqs=0
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            N_seqs+=1
        N_seqs_list.append(str(N_seqs))
    command="python ../classifyFamiliesv5_VennGH_mod_for_candidas_5spp.py "+ str(" ".join(N_seqs_list))+ " %s/synteny_info_output_%s.txt "%(directory,directory) + "%s/%s.fam "%(directory,directory)+"%s/%s.txt "%(directory,directory)+ "%s/%s_FAM_VENN.R "%(directory,directory)+"%s/%s_GENES_VENN.R "%(directory,directory)+"> %s/output_%s_families "%(directory,directory)
        
    print(command)
    subprocess.call(command,shell=True)

make_family_command(species_u, "intergenic")


### run r scripts to make venn diagramms
run_r_scripts = "Rscript intergenic/intergenic_FAM_VENN.R &&\
    Rscript intergenic/intergenic_GENES_VENN.R"
subprocess.call(run_r_scripts,shell=True)    

