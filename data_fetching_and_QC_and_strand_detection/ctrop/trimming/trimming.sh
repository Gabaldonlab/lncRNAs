### NOTE1: before running this step the seqeuncing data must be downloaded
### NOTE2: change the path to the file with adapters
### NOTE3: file_list.txt is a file listing the basename of fastq file, e.g. if a fastq file is named SRRXXXXX_1.fastq.gz, then entry in this file should be SRRXXXXX

while read filename; do
	echo ${filename}
	if [ -e ../${filename}"_2.fastq.gz" ]; then
		java -jar trimmomatic-0.36.jar PE -threads 10 ../${filename}"_1.fastq.gz" ../${filename}"_2.fastq.gz" ${filename}"_tr_1P.fastq.gz" ${filename}"_tr_1U.fastq.gz" ${filename}"_tr_2P.fastq.gz" ${filename}"_tr_2U.fastq.gz" ILLUMINACLIP:path/to/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:49
	else
		java -jar trimmomatic-0.36.jar SE -threads 10 ../${filename}"_1.fastq.gz" ${filename}"_tr_1.fastq.gz" ILLUMINACLIP:path/to/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:49
	fi
done <file_list.txt

mkdir fastqc

fastqc -t 2 *_tr_1.fastq.gz
fastqc -t 2 *_tr_1P.fastq.gz
fastqc -t 2 *_tr_2P.fastq.gz

mv *{zip,html} fastqc
multiqc fastqc/ -o fastqc/
