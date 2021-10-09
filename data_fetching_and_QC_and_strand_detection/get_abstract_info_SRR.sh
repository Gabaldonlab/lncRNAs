while read  -u10 accession;
	do
	experiment_info="$(esearch -db sra -query ${accession} | esummary -format native | xtract -pattern STUDY_ABSTRACT -element STUDY_ABSTRACT)"
	echo ${accession}"!" ${experiment_info}
done 10<all_SRR.txt
