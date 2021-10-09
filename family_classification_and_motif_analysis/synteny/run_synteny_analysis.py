import os
import sys
import subprocess


species_u=["calb_u","ctrop_u","cpar_u","caur_u","cglab_u"]


with open("orthologs.txt","r+") as orthologs: 
	orth = []
	next(orthologs)
	for line in orthologs:
		line=line.rstrip().split("\t")
		orth.append([line[0],line[1],line[2],line[3]])

	
def spp_orth(orth_file):
	orth_dict = {}
	with open(orth_file, "r+") as in_file:
		for line in in_file:
			if not line.startswith("#"):
				line=line.rstrip().split("\t")
				orth_dict[line[0]] = line[3]
	return(orth_dict)
		
calb_caur = spp_orth("C_albicans_SC5314_C_auris_B8441_orthologs.txt")
cpar_caur = spp_orth("C_parapsilosis_CDC317_C_auris_B8441_orthologs.txt")
cglab_caur = spp_orth("C_glabrata_CBS138_C_auris_B8441_orthologs.txt")
			

					
def ctrop_caur_orth(orth_file):
	orth_dict = {}
	with open(orth_file, "r+") as in_file:
		for line in in_file:
			line=line.rstrip().split("\t")
			if line[4] != "---" and line[16] != "---":
				orth_dict[line[4]] = line[16]
	return(orth_dict)		

ctrop_caur = ctrop_caur_orth("cgob_orth.txt")




def generate_five_spp(orth):
	with open("five_spp_orth.txt","w") as five_spp_orth:
		five_spp_orth.write("calb\tctrop\tcpar\tcaur\tcglab\tmyID\n")
		n=1
		for el in orth:
			if el[0] in spp_orth("C_albicans_SC5314_C_auris_B8441_orthologs.txt") and el[1] in ctrop_caur_orth("cgob_orth.txt") and el[2] in spp_orth("C_parapsilosis_CDC317_C_auris_B8441_orthologs.txt") and el[3] in spp_orth("C_glabrata_CBS138_C_auris_B8441_orthologs.txt"):
				five_spp_orth.write("%s\t%s\t%s\t%s\t%s\tmyID%s\n"%(el[0],el[1],el[2],calb_caur[el[0]],el[3],n))
				n+=1
generate_five_spp(orth)

### Run synteny script


def make_family_command(spp_list,directory):
    N_seqs_list=[]
    
    for spp in spp_list:
        if "_u" in spp:
            gene_list_file="gene_lists/%s_protcod_and_u.txt"%(spp.split("_")[0])
        else:
            gene_list_file="gene_lists/%s_protcod_and_u_x.txt"%(spp)
		
        N_seqs = 0
        with open(gene_list_file,"r+") as in_file:
            for line in in_file:
                if "MSTRG" in line:
                    N_seqs+=1
        N_seqs_list.append(str(N_seqs))
        
        
    if "_u" in spp:
		run_synteny = "python ../synteny_nematodesv4GH_mod_for_candidas_5spp.py gene_lists/calb_protcod_and_u.txt gene_lists/ctrop_protcod_and_u.txt gene_lists/cpar_protcod_and_u.txt gene_lists/caur_protcod_and_u.txt gene_lists/cglab_protcod_and_u.txt five_spp_orth.txt %s/synteny_info_output_intergenic.txt 3 3 1 no > %s/output_synteny_intergenic"%(directory,directory)
		classify = "python ../classifyFamiliesv5_VennGH_mod_for_candidas_5spp.py "+ str(" ".join(N_seqs_list))+ " %s/synteny_info_output_intergenic.txt %s/intergenic.fam %s/intergenic.txt %s/intergenic_FAM_VENN.R %s/intergenic_GENES_VENN.R > %s/output_intergenic_families"%(directory,directory,directory,directory,directory,directory)

    
    print(run_synteny)
    print(classify)
    subprocess.call(run_synteny,shell=True)
    subprocess.call(classify,shell=True)
    
    ### run r scripts to make venn diagramms
    run_r_scripts = "Rscript %s/intergenic_FAM_VENN.R && \
    Rscript %s/intergenic_GENES_VENN.R"%(directory,directory)
    subprocess.call(run_r_scripts,shell=True)   

make_family_command(species_u, "intergenic")


