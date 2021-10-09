### See annotations above each command to see its purpose.
### Each command can be executed by uncommenting it at the end of the script.

import sys
import os
import subprocess
from Bio import SeqIO


spp=sys.argv[1]
assemblies = sys.argv[2]




if not os.path.exists("%s/lncRNA_prediction/feelnc"%(spp)):
	os.mkdir("%s/lncRNA_prediction/feelnc"%(spp))
if not os.path.exists("%s/lncRNA_prediction/cpc"%(spp)):
	os.mkdir("%s/lncRNA_prediction/cpc"%(spp))
if not os.path.exists("%s/lncRNA_prediction/counts"%(spp)):
	os.mkdir("%s/lncRNA_prediction/counts"%(spp))
if not os.path.exists("%s/lncRNA_prediction/DE_analysis"%(spp)):
	os.mkdir("%s/lncRNA_prediction/DE_analysis"%(spp))



### extract genes from fasta and gff
extract_genes = "gffread -w %s/lncRNA_prediction/%s_genes.fasta -W -F -g ../reference_genomes_and_annotations/%s.fa \
../reference_genomes_and_annotations/%s.gff"%(spp,spp,spp,spp)

### renamed entries in fasta file
rename_genes = "sed 's/ /|/g' %s/lncRNA_prediction/%s_genes.fasta > %s/lncRNA_prediction/%s_genes_renamed.fasta"%(spp,spp,spp,spp)

### selects only protein coding gene ids
select_prot_cod_ids = "grep -A 1 $'\tgene\t' ../reference_genomes_and_annotations/%s.gff | \
grep $'\tmRNA\t'| cut -f 9| cut -f 2 -d ';'| sed 's/Parent=//g' > %s/lncRNA_prediction/%s_prot_cod_ids.txt"%(spp,spp,spp)

### selects only protein coding sequences
select_prot_cod_fasta = "python ../select_only_prot_cod.py %s/lncRNA_prediction/%s_genes_renamed.fasta %s/lncRNA_prediction/%s_genes_renamed.fasta %s/lncRNA_prediction/%s_only_prot_cod.fasta"%(spp,spp,spp,spp,spp,spp)

### remove prot. coding seqeunces with ambiguous nucleotides
remove_amb_nucl_prot_cod = "python ../remove_seqs_with_amb_nucl.py  \
%s/lncRNA_prediction/%s_only_prot_cod.fasta %s/lncRNA_prediction/feelnc/%s_prot_cod_genes_no_amb_nucl.fasta  1> %s/lncRNA_prediction/feelnc/amb_nuclt_output" %(spp,spp,spp,spp,spp)

### run stringtie merge
merge = "stringtie --merge -o %s/lncRNA_prediction/%s_merged.gtf -g 50 \
                -v -G ../reference_genomes_and_annotations/%s.gff \
                %s"%(spp,spp,spp,assemblies)


### run gffcompare        
compare = "gffcompare -V %s/lncRNA_prediction/%s_merged.gtf -o %s/lncRNA_prediction/%s_merged_compared.gtf \
                -r ../reference_genomes_and_annotations/%s.gff"%(spp, spp, spp, spp,spp)

### select u and x class codes                
non_code = "python ../non_cod.py %s/lncRNA_prediction/%s_merged_compared.gtf.annotated.gtf %s/lncRNA_prediction/unknown_and_antisense_ids.gtf"%(spp,spp,spp)
        

### select only longest transcirpts with cgat
cgat = "source /home/hhovhannisyan/cgat-install/conda-install/bin/activate cgat-s && cgat gtf2gtf \
        --method=filter --filter-method=longest-transcript -I %s/lncRNA_prediction/unknown_and_antisense_ids.gtf \
        > %s/lncRNA_prediction/unknown_and_antisense_ids_longest_transcripts.gtf && conda deactivate"%(spp,spp)


### make fasta file of u and x class codes        
gffread = "gffread -w %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts.fasta -W -F -g ../reference_genomes_and_annotations/%s.fa %s/lncRNA_prediction/unknown_and_antisense_ids_longest_transcripts.gtf"%(spp,spp,spp,spp)
        

### rename 
sed = "sed 's/ /_/g' %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts.fasta > %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed.fasta"%(spp,spp,spp,spp)
     
     
### select transcirpts longer than 200 bp        
select_200 = "python ../select_longer_200.py %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed.fasta \
        %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta"%(spp,spp,spp,spp)
        
### remove sequences with ambigouos nucleotides        
remove_amb_nucl_lncRNA = "python ../remove_seqs_with_amb_nucl.py %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta \
        %s/lncRNA_prediction/feelnc/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta > %s/lncRNA_prediction/removed_trascripts_with_amb_nuclt.txt"%(spp,spp,spp,spp,spp)
        
        
### run FEELnc        
feelnc = "source ~/miniconda2/etc/profile.d/conda.sh && conda activate feelnc &&\
 FEELnc_codpot.pl --outdir='%s/lncRNA_prediction/feelnc/feelnc_codpot_out/' -i %s/lncRNA_prediction/feelnc/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta -a %s/lncRNA_predictionfeelnc/%s_prot_cod_genes_no_amb_nucl.fasta --mode=shuffle \
 && conda deactivate"%(spp,spp,spp,spp,spp)



### select non-coding sequences from FEELnc result
select_non_cod_feelnc= "awk '$11==\"0\"' %s/lncRNA_prediction/feelnc/feelnc_codpot_out/%s_unknown_and_antisense_transcripts_renamed_longer200_no_amb_nucl.fasta_RF.txt| cut -f 1 > %s/lncRNA_prediction/feelnc/feelnc_codpot_out/noncoding_ids_feelnc.txt"%(spp,spp,spp)

### run cpc 
cpc = "~/Softs/cpc-0.9-r2/bin/run_predict.sh %s/lncRNA_prediction/%s_unknown_and_antisense_transcripts_renamed_longer200.fasta %s/lncRNA_prediction/cpc/%s_results.tab %s/lncRNA_prediction/cpc/ %s/lncRNA_prediction/cpc/evidence.txt"%(spp,spp,spp,spp,spp,spp)

### select non-coding sequences from cpc result
select_non_cod_cpc = "awk '$3==\"noncoding\"' %s/lncRNA_prediction/cpc/%s_results.tab | cut -f 1 > %s/lncRNA_prediction/cpc/noncoding_ids_cpc.txt"%(spp,spp,spp)



def find_ovelap():
	with open("%s/lncRNA_prediction/feelnc/feelnc_codpot_out/noncoding_ids_feelnc.txt"%(spp), "r") as feelnc_ids, open("%s/lncRNA_prediction/cpc/noncoding_ids_cpc.txt"%(spp),"r") as cpc_ids, open("%s/lncRNA_prediction/removed_trascripts_with_amb_nuclt.txt"%(spp), "r") as removed_transcipts , open("%s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt"%(spp),"w") as output:
		cpc_noncod=[]
		feelnc_noncod=[]
		cpat_noncod=[]
		
		for line in feelnc_ids:
			line=line.rstrip()
			feelnc_noncod.append(line.upper())
		
		for line in cpc_ids:
			line=line.rstrip()
			cpc_noncod.append(line.upper())
			

		overlap = list(set(feelnc_noncod) & set(cpc_noncod)) 
		for line in overlap:
			output.write("%s\n"%(line))
			
			
		for transcript in removed_transcipts:
			#print transcript
			transcript=transcript.rstrip().upper()
			if transcript in cpc_noncod:
				#print transcript
				output.write("%s\n"%(transcript))


### make a gtf file including coding genes, x and u transcripts. 
helper_gtf = "python ../helper_generate_saf.py %s/lncRNA_prediction/%s_merged_compared.gtf.annotated.gtf %s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt %s/lncRNA_prediction/genes_and_noncod_u_and_x.gtf"%(spp,spp,spp,spp)


### again select only the longest transcripts
cgat2 = "source /home/hhovhannisyan/cgat-install/conda-install/bin/activate cgat-s && cgat gtf2gtf\
        --method=filter --filter-method=longest-transcript -I %s/lncRNA_prediction/genes_and_noncod_u_and_x.gtf  \
        > %s/lncRNA_prediction/genes_and_noncod_u_and_x_longest.gtf  && conda deactivate"%(spp,spp)

### code below generates the saf and bed files using the gtf file obtained above. This whole procedure ensures that all non-coding x and u ids are included in final saf file, and that it does not contain k or p class codes.
generate_saf = "python ../generate_saf_mstrg.py \
%s/lncRNA_prediction/genes_and_noncod_u_and_x_longest.gtf %s/lncRNA_prediction/cpc_feelnc_noncod_ids.txt %s %s/lncRNA_prediction/counts/%s_unsorted.saf %s/lncRNA_prediction/%s_unsorted.bed \
../reference_genomes_and_annotations/%s.gff"%(spp,spp,spp,spp,spp,spp,spp,spp)
sort_saf = "sort -k2,2 -k3,3n %s/lncRNA_prediction/counts/%s_unsorted.saf > %s/lncRNA_prediction/counts/%s.saf"%(spp,spp,spp,spp)
sort_bed = "sort -k1,1 -k2,2n %s/lncRNA_prediction/%s_unsorted.bed > %s/lncRNA_prediction/%s.bed"%(spp,spp,spp,spp)


### generates fasta and bed files for lncRNAs
generate_lncRNA_fasta = "grep 'MSTRG' %s/lncRNA_prediction/%s.bed > %s/lncRNA_prediction/%s_lncRNAs.bed && gffread -w %s/lncRNA_prediction/%s_lncRNAs.fasta -W -F -g ../reference_genomes_and_annotations/%s.fa %s/lncRNA_prediction/%s_lncRNAs.bed"%(spp,spp,spp,spp,spp,spp,spp,spp,spp)

### generate counts for S dataset
bams_for_mapping_bigRNAseq = "ls %s/*subset*bam > %s/lncRNA_predictioncounts/files_for_counting.txt"%(spp,spp)
generate_counts_bigRNAseq = 'while read filename; do prefix_filename=\"$(echo ${filename} | rev | cut -f 2 -d"/" | rev)" && featureCounts -F SAF -p -T 2 -s 2 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts/${prefix_filename}_counts.txt ${filename}; done < %s/lncRNA_prediction/counts/files_for_counting.txt'%(spp,spp,spp,spp)


def generate_counts_public():

	strand_file="../data_fetching_and_QC_and_strand_detection/%s/strand_detection/pseudomapping/strand_info.txt"%(spp)
	with open(strand_file, "r") as strand_file:
		next(strand_file)
		for line in strand_file:
			line=line.rstrip().split("\t")
			if line[3]=="used for transcript reconstruction":
				if "SR" in line[2]:
					featurecounts="featureCounts -F SAF -p -T 2 -s 2 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts/%s_counts.txt ./%s/output_%s/accepted_hits.bam"%(spp,spp,spp,line[0],spp,line[0])
				elif "SF" in line[2]:
					featurecounts="featureCounts -F SAF -p -T 2 -s 1 -a %s/lncRNA_prediction/counts/%s.saf -o %s/lncRNA_prediction/counts/%s_counts.txt ./%s/output_%s/accepted_hits.bam"%(spp,spp,spp,line[0],spp,line[0])
				subprocess.call(featurecounts,shell=True)
				
	
def generate_gc_content_data(output_file, input_file):
	concat_fasta="cat %s/lncRNA_prediction/%s_genes.fasta %s/lncRNA_prediction/%s_lncRNAs.fasta > %s/lncRNA_prediction/%s_genes_and_lncRNAs.fasta"%(spp,spp,spp,spp,spp,spp)
	subprocess.call(concat_fasta, shell=True)
	with open("%s"%(output_file),"w") as gc_output:
		for seq_record in SeqIO.parse("%s"%(input_file), "fasta"):
			name = str(seq_record.id)
			seq = str(seq_record.seq)
			seq=seq.upper()
			gc=float(seq.count("G") + seq.count("C")) / float(len(seq))
			if "|u|" in name:
				class_code="u"
			elif "|x|" in name:
				class_code="x"
			elif "intergenic" in name:
				class_code="inter"
			else:
				class_code="pc"
			
			gc_output.write("%s\t%s\t%s\t%s\n"%(name,class_code,gc,spp))

#### Executing

#subprocess.call(extract_genes, shell=True)
#subprocess.call(rename_genes, shell=True)
#subprocess.call(select_prot_cod_ids, shell=True)
#subprocess.call(select_prot_cod_fasta, shell=True)
#subprocess.call(remove_amb_nucl_prot_cod, shell=True)
#subprocess.call(merge, shell=True)
#subprocess.call(compare, shell=True)
#subprocess.call(non_code, shell=True)
#subprocess.call(cgat, shell=True)
#subprocess.call(gffread, shell=True)
#subprocess.call(sed, shell=True)
#subprocess.call(select_200, shell=True)
#subprocess.call(remove_amb_nucl_lncRNA, shell=True)
#subprocess.call(feelnc, shell=True)
#subprocess.call(select_non_cod_feelnc,shell=True)
#subprocess.call(cpc,shell=True)
#subprocess.call(select_non_cod_cpc,shell=True)
#find_ovelap()
#subprocess.call(helper_gtf, shell=True)
##subprocess.call(cgat2, shell=True)
#subprocess.call(generate_saf,shell=True)
#subprocess.call(sort_saf,shell=True)
#subprocess.call(sort_bed,shell=True)
#subprocess.call(generate_lncRNA_fasta,shell=True)
#subprocess.call(bams_for_mapping_bigRNAseq,shell=True)
#subprocess.call(generate_counts_bigRNAseq,shell=True)
#generate_counts_public()


### this runs gc content for all species
#for spp in ["calb","ctrop","cpar","caur","cglab"]:
#	generate_gc_content_data("%s/lncRNA_prediction/%s_intergenic_gc.tsv"%(spp,spp), \
#	"../%s_intergenic.fa"%(spp))
