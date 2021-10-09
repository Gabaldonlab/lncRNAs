import subprocess
import sys
import os
import itertools
from Bio import SeqIO
import numpy as np



species_u=["calb_u","ctrop_u","cpar_u","caur_u","cglab_u"]
batch_rnafold=open("batch_rnafold.txt","w")


####### Run this part first, generate batch file for rnafold and run it in HPC #######

for iteration in range(1,51):
	#print(iteration)
	if not os.path.exists("iter_"+str(iteration)):
		os.mkdir("iter_"+str(iteration))
	for spp in species:
		
		
		### shuffling

		shuffle_u="fasta-shuffle-letters ../%s_lncRNAs_renamed.fasta iter_%s/%s_lncRNAs_renamed_shuffled.fasta -seed %s"%(spp,iteration,spp,iteration)
		subprocess.call(shuffle_u,shell=True)
		
		### rnafold
		batch_rnafold.write("RNAfold --noPS -i iter_%s/%s_lncRNAs_renamed_shuffled.fasta -j4 > iter_%s/%s_shuffled_lncRNAs_fold.fastb\n"%(iteration,spp,iteration,spp))

		subprocess.call(rnafold_u,shell=True)

batch_rnafold.close()

##### When RNAfold runs in cluster are finished, run the rest ###

### remove extra character from files
remove_brackets='for i in {1..50};do for fold in iter_${i}/{calb,cglab,cpar,caur,ctrop}*fastb;do fbname=$(basename "${fold}") && cut -f1 -d " " ${fold} > iter_${i}/final_${fbname};done;done'
#subprocess.call(remove_brackets,shell=True)



### parse RNAfold files to fix Beagle bug
def fix_RNAfold(species):
	for iteration in range(1,51):
		for spp in species:
			for spp_file in [spp,spp+"_u"]:
				spp_dict={}
				with open("iter_%s/final_%s_shuffled_lncRNAs_fold.fastb"%(iteration,spp_file),"r+") as fold_file, open("iter_%s/final_%s_shuffled_lncRNAs_fold_for_beagle.fastb"%(iteration,spp_file),"w") as output:
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

#fix_RNAfold(species)



#### this will generate batch file for beagle, run it in HPC

def generate_beagle_batch(spp):
	if "_u" in spp[1]:
		batch_beagle=open("batch_beagle_u.txt","w")
		
	for iteration in range(1,51):
		all_pairwise_comparisons_without_repeats_and_without_doubles=list(itertools.combinations(spp, 2))
		print(all_pairwise_comparisons_without_repeats_and_without_doubles)


		for comparison in all_pairwise_comparisons_without_repeats_and_without_doubles:
			
			spp1=comparison[0]
			spp2=comparison[1]
			
			## make directory for each comparison
			mkdir="mkdir -p iter_%s/%s_to_%s/"%(iteration,spp1,spp2)
			subprocess.call(mkdir, shell=True)
			
			## splitting the file for the spp1
			split_file = "split -l 3 iter_%s/final_%s_shuffled_lncRNAs_fold_for_beagle.fastb -d -a4 %s_split_ --additional-suffix=.fastb && mv %s_split* iter_%s/%s_to_%s/"%(iteration,spp1,spp1,spp1,iteration,spp1,spp2)
			subprocess.call(split_file, shell=True)
			for root, dirs, files in os.walk("iter_%s/%s_to_%s/"%(iteration,spp1,spp2)):
				for f in files:
					if "split" in f:
						print("java -jar -Xloggc:%s_to_%s.log ../../Beagle_v0.2.jar -input1 ../iter_%s/%s_to_%s/%s -input2 ../iter_%s/final_%s_shuffled_lncRNAs_fold_for_beagle.fastb -c 2 -l true -pValue true -outfile ../iter_%s/%s_to_%s/%s_aligned.txt\n"%(spp1,spp2,iteration,spp1,spp2,f,iteration,spp2,iteration,spp1,spp2,f))
						batch_beagle.write("java -jar -Xloggc:%s_to_%s.log ../../Beagle_v0.2.jar -input1 ../iter_%s/%s_to_%s/%s -input2 ../iter_%s/final_%s_shuffled_lncRNAs_fold_for_beagle.fastb -c 2 -l true -pValue true -outfile  ../iter_%s/%s_to_%s/%s_aligned.txt\n"%(spp1,spp2,iteration,spp1,spp2,f,iteration,spp2,iteration,spp1,spp2,f))
	batch_beagle.close()

#generate_beagle_batch(species_u)

### The following lines are specific for the HPC at Barcelona Supercomputing Centre. Adjust to your needs.


def prepare_for_MN_launch(spp):
	if "_u" in spp[1]:
		if not os.path.exists("launch_beagle_u"):
			os.mkdir("launch_beagle_u")
		split_u='cp batch_beagle_u.txt launch_beagle_u/ && cd launch_beagle_u/ && split -l 15000 batch_beagle_u.txt -d -a5 job_ && for j in job_*; do echo "/apps/GREASY/2.2/INTEL/IMPI/bin/greasy ${j}"; done > greasy_array_beagle.txt && cd ..'
		subprocess.call(split_u,shell=True)


#prepare_for_MN_launch(species_u)


### remove failed files
def rename_files_to_exclude(folder):
	for root, dirs, files in os.walk(folder):
		for f in files:
			if "rst" in f:
				with open(folder+f,"r") as  rst:
					for l in rst:
						if l.startswith("java"):
							#print("%s\t%s\n"%(f,l))
							file_to_exclude=l.split(" ")[-1].replace("../","").rstrip()
							rename="mv %s %s_EXCLUDED"%(file_to_exclude,file_to_exclude)
							print(rename)
							subprocess.call(rename,shell=True)

#rename_files_to_exclude("launch_beagle_u/")


def generate_classification(spp):
	percentage_of_tr_in_fam=[]
	n_iter=0
	n_random_more_than_real=0
	for iteration in range(1,51):
		n_iter+=1
		#print(iteration)
		all_pairwise_comparisons_without_repeats_and_without_doubles=list(itertools.combinations(spp, 2))
		#print(all_pairwise_comparisons_without_repeats_and_without_doubles)


		for comparison in all_pairwise_comparisons_without_repeats_and_without_doubles:
			
			spp1=comparison[0]
			spp2=comparison[1]
			
			## concatenate files for each comparison
			make_header_files='cat iter_%s/%s_to_%s/*aligned.txt | grep "^>MSTR" > iter_%s/%s_to_%s/%s_to_%s_headers.txt'%(iteration,spp1,spp2,iteration,spp1,spp2,spp1,spp2)
			#subprocess.call(make_header_files, shell=True)
			
		if "u" in spp[1]:
			concat_headers='cat iter_%s/{calb,cpar,ctrop,caur,cglab}_u_to*/*headers.txt > iter_%s/all_headers_u.txt'%(iteration,iteration)

		#subprocess.call(concat_headers, shell=True)
		
		
		if "u" in spp[1]:
			file_for_cinta_script = "python ../parse_beagle.py iter_%s/all_headers_u.txt intergenic iter_%s"%(iteration,iteration)

		#subprocess.call(file_for_cinta_script,shell=True)
		
		if "u" in spp[1]:
			classify='sed "s/_shuf//g" iter_%s/results_for_cintas_script_intergenic.txt > iter_%s/results_for_cintas_script_intergenic_final.txt && python ../../classifyFamiliesv5_VennGH_mod_for_candidas_5spp.py 1459 1581 1499 842 449 iter_%s/results_for_cintas_script_intergenic_final.txt iter_%s/intergenic.fam iter_%s/intergenic.txt iter_%s/intergenic_FAM_VENN.R iter_%s/intergenic_GENES_VENN.R >iter_%s/output_intergenic_families'%(iteration,iteration,iteration,iteration,iteration,iteration,iteration,iteration)

		#subprocess.call(classify, shell=True)
		
		real=5724
		
		if "u" in spp[1]:
			fam = subprocess.check_output('wc -l iter_%s/intergenic.fam'%(iteration), shell = True)
			N_trans_in_fam = float(fam.split(" ")[0])
			#print(N_trans_in_fam/float(5830)*100)
			if N_trans_in_fam >= real:
				n_random_more_than_real+=1
				#print(iteration)
			percentage_of_tr_in_fam.append(N_trans_in_fam/float(5830)*100)	
			print(iteration, N_trans_in_fam/float(5830)*100)

		else:
			fam = subprocess.check_output('wc -l iter_%s/all.fam'%(iteration), shell = True)
			N_trans_in_fam = float(fam.split(" ")[0])
			print(N_trans_in_fam/float(10741)*100)
			percentage_of_tr_in_fam.append(N_trans_in_fam/float(10741)*100)	
	print(percentage_of_tr_in_fam)
	with open("reshuffle_final_result.txt","w") as final_result:
		print("%s simulation sets are more/equal than real %s%%, pvalue = %s"%(n_random_more_than_real,round(float(real)/float(5830),3)*100,float(n_random_more_than_real)/float(n_iter)))
		print("On average %s%% (+- %s)  of randomly reshufled lncRNAs form families"%(round(np.mean(percentage_of_tr_in_fam),2),round(np.std(percentage_of_tr_in_fam),2)))

		final_result.write("%s simulation sets are more/equal than real %s%%, pvalue = %s\n"%(n_random_more_than_real,round(float(real)/float(5830),3),float(n_random_more_than_real)/float(n_iter)))
		final_result.write("On average %s%% (+- %s)  of randomly reshufled lncRNAs form families"%(round(np.mean(percentage_of_tr_in_fam),2),round(np.std(percentage_of_tr_in_fam),2)))

#generate_classification(species_u)
