spp_list=["calb","cglab","cpar","caur","ctrop"]


for spp in spp_list:
	saf_file="/home/hhovhannisyan/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/%s/counts/%s.saf"%(spp,spp)
	genelist_protcod_u="/home/hhovhannisyan/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/synteny_clustering/gene_lists/%s_protcod_and_u.txt"%(spp)
	genelist_protcod_u_x="/home/hhovhannisyan/users/tg/current/hhovhannisyan/jena_project/lncRNAs_bigRNAseq_and_public/synteny_clustering/gene_lists/%s_protcod_and_u_x.txt"%(spp)

	with open(saf_file, "r") as saf, open(genelist_protcod_u, "w") as genes_and_u, open(genelist_protcod_u_x,"w") as genes_and_u_x:
		#next(saf)
		for line in saf:
			line=line.rstrip().split("\t")
			if "|x|" in line[0]:
				genes_and_u_x.write("%s\n"%(line[0]))
			elif "|u|" in line[0]:
				genes_and_u.write("%s\n"%(line[0]))
				genes_and_u_x.write("%s\n"%(line[0]))
			else:
				genes_and_u.write("%s\n"%(line[0].split("|")[0]))
				genes_and_u_x.write("%s\n"%(line[0].split("|")[0]))
