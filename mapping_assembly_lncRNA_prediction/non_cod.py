import sys

with open(sys.argv[1], "r+") as gtf, open(sys.argv[2],"w") as output:
	for line in gtf:
		line=line.rstrip().split("\t")
		if line[2] == "transcript":
			#print line[-1]
			if ('class_code "x"' in line[-1]) or ('class_code "u"' in line[-1]):
				var="+"
				output.write("%s\n"%("\t".join(line)))
				#print line[-1]
			else:
				var="-"
		elif line[2]=="exon":
			if var=="+":
				#print "\t".join(line)
				output.write("%s\n"%("\t".join(line)))
	
