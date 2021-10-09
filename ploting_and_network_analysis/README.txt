This directory contains scripts and data files to perform analysis and plotting of most of the results.
The main script is called lncRNA_results.R. The libraries and versions necessary for the analysis are listed in sessionInfo.txt. We recommend to run this script step by step in Rstudio.
Each directory contains all  necessary files for the analysis per species.
In particular:

*_chrNameLength.txt - length and names of chromosomes
*_gc.tsv  - GC content data
*_lncRNAs.bed - bed file of lncRNAS
infection_lncRNAs_*_24_24c.txt - list of infection-specific lncRNAS
*_iprscan.out - file for performing PFAM enrichment
*_distance.bed - distance to the telomeres data
*_go.txt - GO tables
*_prot_cod_ids.txt - protein coding ids
saturation_data_*.txt - data of saturation plot analysis
*.bed - bed ifle of lncRNAs and other features
*_expression.txt  - read count data
*_intergenic_gc.tsv - GC content of intergenic regions
*_SraRunTable.txt - metadata from SRA
SRR_SRP_correspondence_*.txt - correspondece between SRR and SRP accession numers

caur_strain_orthologs.tsv - orthogs of C. auris strains
{caur,ctrop}_to_discard.txt - lncRNAs of C. auris or C. tropicalis to discard (potential tRNAs)
