### NOTE1: file_list.txt is a file listing the basename of trimmed fastq files, e.g. if a fastq file is named SRRXXXXX_tr_1P.fastq.gz, then entry in this file should be SRRXXXXX

while read filename;
do
NAME=$(basename ${filename})
echo "Started mapping $filename"
if [ -e ../../trimming/${filename}"_tr_2P.fastq.gz" ]; then
	salmon quant -i  ../ref_gen/calb -l A \
	-1 ../../trimming/${filename}_tr_1P.fastq.gz -2 ../../trimming/${filename}_tr_2P.fastq.gz \
	-g ../ref_gen/calb.gff \
	-p 11 -o ${NAME} --gcBias
else
    salmon quant -i  ../ref_gen/calb -l A \
    -r ../../trimming/${filename}_tr_1.fastq.gz \
    -g ../ref_gen/calb.gff \
    -p 11 -o ${NAME} --gcBias
fi
done <file_list.txt

