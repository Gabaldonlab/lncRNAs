import subprocess
import os
import glob






def run_bedtools():
    stringtie_trinity_uniq_prop = open("./unique_transcirpts_trinity_vs_stringtie.txt","w")
    for spp in ["calb","ctrop","cpar","caur","cglab"]:
        for sample_path in glob.glob("./%s/*_trinity"%(spp)):
            sample = os.path.basename(sample_path).replace("_trinity","")

            generate_trinity_genes_gff = """awk '$3=="gene"' %s/%s_trinity.gff > %s/%s_trinity_genes.gff"""%(sample_path,sample,sample_path,sample)
            #print(generate_trinity_genes_gff)
            #subprocess.call(generate_trinity_genes_gff,shell=True)

            generate_stringtie_genes_gff = """awk '$3=="gene"' %s/%s_stringtie.gff > %s/%s_stringtie_genes.gff"""%(sample_path,sample,sample_path,sample)
            #print(generate_stringtie_genes_gff)
            #subprocess.call(generate_stringtie_genes_gff,shell=True)



            ### Trinity cmds
            only_trinity_cmd = "bedtools intersect -v -a %s/%s_trinity_genes.gff -b %s/%s_stringtie_genes.gff | wc -l"%(sample_path,sample,sample_path,sample)
            only_trinity = subprocess.check_output(only_trinity_cmd,shell=True)

            total_trinity_cmd = "cat %s/%s_trinity_genes.gff | wc -l "%(sample_path,sample)
            total_trinity = subprocess.check_output(total_trinity_cmd,shell=True)

            print("%s\t%s\t%s\t%s\n"%(int(only_trinity),int(total_trinity), "Trinity",spp))
            stringtie_trinity_uniq_prop.write("%s\t%s\t%s\t%s\n"%(int(only_trinity),int(total_trinity), "Trinity",spp))



            stringtie_in_trinity = "bedtools intersect -a %s/%s_trinity_genes.gff  -b %s/%s_stringtie_genes.gff -c > %s/%s_stringtie_in_trinity.txt"%(sample_path,sample,sample_path,sample,sample_path,sample)
            subprocess.call(stringtie_in_trinity,shell=True)


            ##################
            ### Stringtie cmds
            only_stringtie_cmd = "bedtools intersect -v -a %s/%s_stringtie_genes.gff -b %s/%s_trinity_genes.gff | wc -l"%(sample_path,sample,sample_path,sample)
            only_stringtie = subprocess.check_output(only_stringtie_cmd,shell=True)

            total_stringtie_cmd = "cat %s/%s_stringtie_genes.gff | wc -l "%(sample_path,sample)
            total_stringtie = subprocess.check_output(total_stringtie_cmd,shell=True)

            print("%s\t%s\t%s\t%s\n"%(int(only_stringtie),int(total_stringtie), "Stringtie",spp))
            stringtie_trinity_uniq_prop.write("%s\t%s\t%s\t%s\n"%(int(only_stringtie),int(total_stringtie), "Stringtie",spp))

            trinity_in_stringtie = "bedtools intersect -a %s/%s_stringtie_genes.gff  -b %s/%s_trinity_genes.gff -c > %s/%s_trinity_in_stringtie.txt"%(sample_path,sample,sample_path,sample,sample_path,sample)
            subprocess.call(trinity_in_stringtie,shell=True)

        merge_trinity_in_stringtie = " cat %s/*trinity/*_trinity_in*txt > %s/all_trinity_in_stringtie.txt"%(spp,spp)
        subprocess.call(merge_trinity_in_stringtie,shell=True)

        merge_stringtie_in_trinity = " cat %s/*trinity/*_stringtie_in*txt > %s/all_stringtie_in_trinity.txt"%(spp,spp)
        subprocess.call(merge_stringtie_in_trinity,shell=True)

    stringtie_trinity_uniq_prop.close()


run_bedtools()



def run_bedtools_ux():
    stringtie_trinity_uniq_prop = open("./unique_transcirpts_trinity_vs_stringtie_ux.txt","w")
    for spp in ["calb","ctrop","cpar","caur","cglab"]:
        for sample_path in glob.glob("./%s/*_trinity"%(spp)):
            sample = os.path.basename(sample_path).replace("_trinity","")
            print(sample_path)
            ######
            generate_trinity_transcripts_gff_u = """awk '$3=="transcript"' %s/%s_trinity_compare.annotated.gtf | grep -e 'class_code "u"' > %s/%s_trinity_compare.annotated_transcript_u.gtf"""%(sample_path,sample,sample_path,sample)
            generate_trinity_transcripts_gff_x = """awk '$3=="transcript"' %s/%s_trinity_compare.annotated.gtf | grep -e 'class_code "x"' > %s/%s_trinity_compare.annotated_transcript_x.gtf"""%(sample_path,sample,sample_path,sample)

            subprocess.call(generate_trinity_transcripts_gff_u,shell=True)
            subprocess.call(generate_trinity_transcripts_gff_x,shell=True)
            ######


            generate_stringtie_transcripts_gff_u = """awk '$3=="transcript"' %s/%s_stringtie_compare.annotated.gtf | grep -e 'class_code "u"' > %s/%s_stringtie_compare.annotated_transcripts_u.gtf"""%(sample_path,sample,sample_path,sample)
            generate_stringtie_transcripts_gff_x = """awk '$3=="transcript"' %s/%s_stringtie_compare.annotated.gtf | grep -e 'class_code "x"' > %s/%s_stringtie_compare.annotated_transcripts_x.gtf"""%(sample_path,sample,sample_path,sample)

            subprocess.call(generate_stringtie_transcripts_gff_u,shell=True)
            subprocess.call(generate_stringtie_transcripts_gff_x,shell=True)


            ### Trinity cmds
            for tr in ["u","x"]:
                only_trinity_cmd = "bedtools intersect -v -a %s/%s_trinity_compare.annotated_transcript_%s.gtf -b %s/%s_stringtie_compare.annotated_transcripts_%s.gtf | wc -l"%(sample_path,sample,tr,sample_path,sample,tr)
                only_trinity = subprocess.check_output(only_trinity_cmd,shell=True)

                total_trinity_cmd = "cat %s/%s_trinity_compare.annotated_transcript_%s.gtf | wc -l "%(sample_path,sample,tr)
                total_trinity = subprocess.check_output(total_trinity_cmd,shell=True)

                #print("%s\t%s\t%s\t%s\n"%(int(only_trinity),int(total_trinity), "Trinity",spp))
                stringtie_trinity_uniq_prop.write("%s\t%s\t%s\t%s\t%s\n"%(int(only_trinity),int(total_trinity), "Trinity",spp,tr))

                stringtie_in_trinity = "bedtools intersect -a %s/%s_trinity_compare.annotated_transcript_%s.gtf  -b %s/%s_stringtie_compare.annotated_transcripts_%s.gtf -c > %s/%s_stringtie_in_trinity_%s.txt"%(sample_path,sample,tr,sample_path,sample,tr,sample_path,sample,tr)
                subprocess.call(stringtie_in_trinity,shell=True)

                ### Stringtie cmds
                only_stringtie_cmd = "bedtools intersect -v -a %s/%s_stringtie_compare.annotated_transcripts_%s.gtf -b %s/%s_trinity_compare.annotated_transcript_%s.gtf  | wc -l"%(sample_path,sample,tr,sample_path,sample,tr)
                only_stringtie = subprocess.check_output(only_stringtie_cmd,shell=True)

                total_stringtie_cmd = "cat %s/%s_stringtie_compare.annotated_transcripts_%s.gtf | wc -l "%(sample_path,sample,tr)
                total_stringtie = subprocess.check_output(total_stringtie_cmd,shell=True)

                stringtie_trinity_uniq_prop.write("%s\t%s\t%s\t%s\t%s\n"%(int(only_stringtie),int(total_stringtie),"Stringtie",spp,tr))

                trinity_in_stringtie = "bedtools intersect -a %s/%s_stringtie_compare.annotated_transcripts_%s.gtf  -b %s/%s_trinity_compare.annotated_transcript_%s.gtf -c > %s/%s_trinity_in_stringtie_%s.txt"%(sample_path,sample,tr,sample_path,sample,tr,sample_path,sample,tr)
                subprocess.call(trinity_in_stringtie,shell=True)

         merge_trinity_in_stringtie = " cat %s/*trinity/*_trinity_in*txt > %s/all_trinity_in_stringtie.txt"%(spp,spp)
         subprocess.call(merge_trinity_in_stringtie,shell=True)


    stringtie_trinity_uniq_prop.close()


run_bedtools_ux()
