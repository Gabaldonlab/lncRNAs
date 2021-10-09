For each species, this folder contains scripts for fetching the sra data, converting it to fastq format, performing quality control, trimming, and detecting strand information.

Required software:
sra-tools v. 2.8.0 or higher
fastqc v. 0.11.8 or other
multiqc v. 1.0 or other
trimmomatic v. 0.36
RSEM v. 1.3 or higher
salmon v. 0.8.1 or higher
python 2 or 3


To perform the analysis for each species:
1. run prefetch scripts
2. run trimming scripts
3. in strand_detection run ref_gen/reference.sh and ref_gen/salmon_index.sh, then run pseudomapping/pseudomapping.sh and pseudomapping/generate_strand_info.py


Additionally, get_abstract_info_SRR.sh script will fetch the experimental data for all analyzed samples.
