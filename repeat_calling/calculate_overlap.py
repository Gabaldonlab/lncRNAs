import sys
import subprocess
import os


def make_repeat_bed(species):
    for spp in species:

        with open("./%s/reference_genome_dir/reference_genome.fasta.repeats.tab"%(spp), "r+") as rep_modeler, open("./%s/reference_genome_dir/%s_repeats.bed"%(spp,spp),"w") as bed:
            next(rep_modeler)
            for line in rep_modeler:
                line = line.rstrip().split()
                #print("%s\t%s\t%s\n"%(line[3],line[2],line[4]))
                if "LINE" in line[14] or "Low_complexity" in line[14] or "LTR" in line[14] or "Simple_repeat" in line[14] or "SINE" in line[14] or "Unknown" in line[14]:
                    bed.write("%s\t%s\t%s\t%s|%s\n"%(line[3],line[2],line[4],line[0],line[14]))

#make_repeat_bed(["calb","ctrop","cpar","caur","cglab"])


def run_bedtools(species):
    with open("repeat_overlap_stats.txt","w") as stats:
        stats.write("%s\t%s\t%s\t%s\t%s\n"%("Species","Type","N","Total","Class_code"))
        for spp in species:
			if spp == "ctrop":
				bedtools = """bedtools intersect -wao -a ../../%s/%s_lncRNAs_notrnas.bed -b ./%s/reference_genome_dir/%s_repeats.bed -F 0.5 > ./%s/%s_lncRNA_repeat_intersection.txt """%(spp,spp,spp,spp,spp,spp)
			elif spp == "caur":
				bedtools = """bedtools intersect -wao -a ../../%s/%s_lncRNAs_notrnas.bed -b ./%s/reference_genome_dir/%s_repeats.bed -F 0.5 > ./%s/%s_lncRNA_repeat_intersection.txt """%(spp,spp,spp,spp,spp,spp)
			else:
				bedtools = """bedtools intersect -wao -a ../../%s/%s_lncRNAs.bed -b ./%s/reference_genome_dir/%s_repeats.bed -F 0.5 > ./%s/%s_lncRNA_repeat_intersection.txt """%(spp,spp,spp,spp,spp,spp)
			print(bedtools)
			subprocess.call(bedtools,shell=True)

			bedtools_pc = """bedtools intersect -wao -a ../../%s/%s_prot_cod.bed -b ./%s/reference_genome_dir/%s_repeats.bed -F 0.5 > ./%s/%s_pc_repeat_intersection.txt """%(spp,spp,spp,spp,spp,spp)
			print(bedtools_pc)
			subprocess.call(bedtools_pc,shell=True)


			test_dict = {}
			with open("./%s/%s_lncRNA_repeat_intersection_selected_repeats.txt"%(spp,spp)) as overlap:
				for line in overlap:
					line=line.rstrip().split("\t")
					if not line[3] in test_dict:
						test_dict[line[3]] = [line[8]]
					else:
						test_dict[line[3]].append(line[8])

			N_multirepeat = 0
			for k,v in test_dict.items():
				if len(v)>1:
					multi_repeat = "|".join(list(set([x.split("|")[1] for x in v])))
                    #print(spp)
					if len(multi_repeat.split("|")) > 1:
						#print(multi_repeat)
						if len(multi_repeat.split("|")) == 2:
							N_multirepeat+=1
						elif len(multi_repeat.split("|")) == 3:
							N_multirepeat+=2
						elif len(multi_repeat.split("|")) == 4:
							N_multirepeat+=3
						else:
							print("CRAZY FUCK",k,v)


			for class_code in ["u","x"]:

				Total = subprocess.check_output("""grep "|%s|" ../../%s/%s_lncRNAs.bed  | wc -l"""%(class_code,spp,spp), shell=True)
                #Total_i = subprocess.check_output("""grep "|u|" ../../%s/%s_lncRNAs.bed  | wc -l"""%(spp,spp), shell=True)

                ### different repeats
				N_lncRNAs_overlaping_LINE=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep -e "LINE" -e "SINE"| cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
				stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"LINE/SINE",N_lncRNAs_overlaping_LINE.rstrip(),Total.rstrip(),class_code))

				N_lncRNAs_overlaping_Low_complex=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "Low_complexity"| cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
				stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"Low complexity repeat",N_lncRNAs_overlaping_Low_complex.rstrip(),Total.rstrip(),class_code))

				N_lncRNAs_overlaping_LTR=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "LTR"| cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
				stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"LTR",N_lncRNAs_overlaping_LTR.rstrip(),Total.rstrip(),class_code))

				N_lncRNAs_overlaping_Simple_repeat=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "Simple_repeat"| cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
				stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"Simple repeat",N_lncRNAs_overlaping_Simple_repeat.rstrip(),Total.rstrip(),class_code))

                #N_lncRNAs_overlaping_SINE=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "SINE"| cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
                #stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"SINE",N_lncRNAs_overlaping_SINE.rstrip(),Total.rstrip(),class_code))

				N_lncRNAs_overlaping_Unknown=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "Unknown"| cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
				stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"Unknown",N_lncRNAs_overlaping_Unknown.rstrip(),Total.rstrip(),class_code))


				N_lncRNAs_not_overlaping_repeats=subprocess.check_output("""awk '$10==0' ./%s/%s_lncRNA_repeat_intersection.txt | cut -f 4|sort| uniq| grep "|%s|" | wc -l"""%(spp,spp,class_code), shell=True)
				stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"No overlap",N_lncRNAs_not_overlaping_repeats.rstrip(),Total.rstrip(),class_code))


            ### for protein coding genes

			Total = subprocess.check_output("""cat ../../%s/%s_prot_cod.bed  | wc -l"""%(spp,spp), shell=True)
            #Total_i = subprocess.check_output("""grep "|u|" ../../%s/%s_lncRNAs.bed  | wc -l"""%(spp,spp), shell=True)

            ### different repeats
			N_pc_overlaping_LINE=subprocess.check_output("""awk '$10!=0' ./%s/%s_pc_repeat_intersection.txt | grep -e "LINE" -e "SINE"| cut -f 4|sort| uniq| wc -l"""%(spp,spp), shell=True)
			stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"LINE/SINE",N_pc_overlaping_LINE.rstrip(),Total.rstrip(),"pc"))

			N_pc_overlaping_Low_complex=subprocess.check_output("""awk '$10!=0' ./%s/%s_pc_repeat_intersection.txt | grep "Low_complexity"| cut -f 4|sort| uniq | wc -l"""%(spp,spp), shell=True)
			stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"Low complexity repeat",N_pc_overlaping_Low_complex.rstrip(),Total.rstrip(),"pc"))

			N_pc_overlaping_LTR=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "LTR"| cut -f 4|sort| uniq| wc -l"""%(spp,spp), shell=True)
			stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"LTR",N_pc_overlaping_LTR.rstrip(),Total.rstrip(),"pc"))

			N_pc_overlaping_Simple_repeat=subprocess.check_output("""awk '$10!=0' ./%s/%s_pc_repeat_intersection.txt | grep "Simple_repeat"| cut -f 4|sort| uniq|  wc -l"""%(spp,spp), shell=True)
			stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"Simple repeat",N_pc_overlaping_Simple_repeat.rstrip(),Total.rstrip(),"pc"))

            #N_pc_overlaping_SINE=subprocess.check_output("""awk '$10!=0' ./%s/%s_pc_repeat_intersection.txt | grep "SINE"| cut -f 4|sort| uniq| wc -l"""%(spp,spp), shell=True)
            #stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"SINE",N_pc_overlaping_SINE.rstrip(),Total.rstrip(),"pc"))

			N_pc_overlaping_Unknown=subprocess.check_output("""awk '$10!=0' ./%s/%s_lncRNA_repeat_intersection.txt | grep "Unknown"| cut -f 4|sort| uniq| wc -l"""%(spp,spp), shell=True)
			stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"Unknown",N_pc_overlaping_Unknown.rstrip(),Total.rstrip(),"pc"))


			N_pc_not_overlaping_repeats=subprocess.check_output("""awk '$10==0' ./%s/%s_pc_repeat_intersection.txt | cut -f 4|sort| uniq| wc -l"""%(spp,spp), shell=True)
			stats.write("%s\t%s\t%s\t%s\t%s\n"%(spp,"No overlap",N_pc_not_overlaping_repeats.rstrip(),Total.rstrip(),"pc"))




			print(spp, N_multirepeat)
run_bedtools(["calb","ctrop","cpar","caur","cglab"])
