while read file; do
    	echo "downloading $file"
    	prefetch -v $file
		mv ~/ncbi/public/sra/${file}.sra ./
		fastq-dump --split-files ${file}.sra
		rm ${file}.sra
		gzip *.fastq
done < cglab_accessions.txt

mkdir ./fastqc
fastqc *.fastq.gz -t 10
mv *.fastq.gz ./fastqc
mv *{zip,html} ./fastqc
multiqc ./fastqc/ -o ./fastqc/
