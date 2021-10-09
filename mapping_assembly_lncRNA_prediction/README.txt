For each species, this folder contains scripts for mapping data to reference genome, asssembling transcriptome, making a unified transcriptome, predict lncRNAs, calculate expression levels and GC content 

Required software:

tophat2 v. 2.1.1 or higher
stringtie v. 1.3.3b or higher
samtools v. 1.9 or higher
gffread v. 0.11.4 or higher
gffcompare v. 0.11.2 or higher 
cgat v 0.3.2 or higher
FEELnc v. 0.1.1 or higher
cpc v. 0.9 or higher
featurecounts v. 1.6.4 or higher
python2
Biopython
Trinity v. 2.8.5
GMAP v. 2016-11-07


To perform the analysis for each species:
1. Run S_dataset_commands.txt and B_dataset_commands.txt 
2. When finished, run run_lncRNA_predictions.sh. See remarks and comments in the script.
3. (Optional) perform comparisons with Trinity assemblies. In the folder ./trinity follow the README.txt.
