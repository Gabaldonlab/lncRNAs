# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 10:55:13 2016

@author: cpegueroles
"""
#usage: python classifyFamiliesv5_VennGH.py calb ctrop cpar cglab 4spv4.out 4spv4_Families.fam 4spv4_Families.txt 4spv4_FAMvenn.R 4spv4_GENESvenn.R >4spv4_Families.counts
#ex: 4spv4.out
#calb    cpar	XLOC_014828	XLOC_024666
#cpar   calb	XLOC_024666	XLOC_014828


import sys,os,re

def renameGenes(in1):#ex: 4spv4_6cluster3min.out
    #renames genes adding the species; ex: XLOC_010119 -> XLOC_010119ctrop
    #creates a dictionary with sp1 (key) pointing to sp2 (val)
    dictTemp={}
    for row in open(in1, 'r').readlines():
        row = row.rstrip().split('\t')
        key= row[2]#+row[0]
        val=[]
        if not key in dictTemp:
            dictTemp[key]=val
        dictTemp[key].append(row[2]) #+row[0])
        dictTemp[key].append(row[3]) #+row[1])
    #print(dictTemp)
    return dictTemp
    
##################################################
calb= sys.argv[1] #number of expressed (protein coding/lncRNA) genes
ctrop= sys.argv[2] #number of expressed (protein coding/lncRNA) genes
cpar= sys.argv[3] #number of expressed (protein coding/lncRNA) genes
caur=sys.argv[4]
cglab= sys.argv[5] #number of expressed (protein coding/lncRNA) genes
in1= sys.argv[6] #out file from synteny_nematodesv4GH.py; ex: 4spv4_6cluster3min.out
out= open(sys.argv[7], 'w') #.fam
out2= open(sys.argv[8], 'w') #.txt
out3= open(sys.argv[9], 'w') #number of families
out4= open(sys.argv[10], 'w') #number of genes

#to group into families
dictTemp1= renameGenes(in1)

dictTemp2= dictTemp1
#print(dictTemp1)
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1

dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1
           
mylist=[val for val in dictTemp1.values()]
mylistSorted=[]
for x in mylist:
    mylistSorted.append(sorted(x)) 
uniq_mylist = [list(t) for t in set(sorted(map(tuple, mylistSorted)))]
     
#print(uniq_mylist)
dictFam={}
#dictFam structure: fam -> {sp1:[,], sp2:[,], sp3:[,], sp4:[,]}
i=1
for x in uniq_mylist:
    key=i
    val=x
    dictFam[key]=val
    i=i+1

calbm=0; ctropm=0; cparm=0; caurm=0; cglabm=0; #number of classified genes
for key, val in dictFam.iteritems():
    for x in val:
        out.write("fam%s\t%s\n" % (key, x))
        if re.search("calb", x):
            calbm=calbm+1
        if re.search("ctrop", x):
            ctropm=ctropm+1
        if re.search("cpar", x):
            cparm=cparm+1
        if re.search("cglab", x):
            cglabm=cglabm+1
        if re.search("caur", x):
			caurm=caurm+1

#to count the number of species for gene family        
albpar=0; albtrop=0;albaur=0; albglab=0; partrop=0; paraur=0; parglab=0; tropaur=0; tropglab=0; aurglab=0
albpartrop=0; albparaur=0; albparglab=0; albtropaur=0; albtropglab=0;albaurglab=0;partropaur=0; partropglab=0; paraurglab=0;tropaurglab=0
albpartropaur=0; albpartropglab=0; albparaurglab=0; albtropaurglab=0; partropaurglab=0
albpartropaurglab=0

genesalbpar=0; genesalbtrop=0; genesalbaur=0; genesalbglab=0; genespartrop=0; genesparaur=0; genesparglab=0; genestropaur=0; genestropglab=0; genesaurglab=0
genesalbpartrop=0; genesalbparaur=0; genesalbparglab=0; genesalbtropaur=0; genesalbtropglab=0; genesalbaurglab=0; genespartropaur=0; genespartropglab=0; genesparaurglab=0; genestropaurglab=0
genesalbpartropaur=0; genesalbpartropglab=0; genesalbparaurglab=0; genesalbtropaurglab=0; genespartropaurglab=0
genesalbpartropaurglab=0

#genesalbpar=0; genesalbtrop=0; genesalbglab=0; genespartrop=0; genesparglab=0; genestropglab=0
#genesalbpartrop=0; genesalbparglab=0; genesalbtropglab=0; genespartropglab=0
#genesalbpartropglab=0

albF=0; tropF=0; parF=0; aurF=0; glabF=0
albG=0; tropG=0; parG=0; aurG=0; glabG=0
for key, val in dictFam.iteritems():
    print(key,val)
    mySp=[]
    for x in val:
        print(x)
        x=x.split('|c')
        sp = x[1]
        mySp.append(sp)
    out2.write("fam%s\t%s\n" % (key,'\t'.join(list(set(mySp)))))
    # 5 spp
    if 'alb' in mySp and 'trop' in mySp and 'par' in mySp and 'aur' in mySp and 'glab' in mySp:
        albpartropaurglab +=1
        genesalbpartropaurglab = genesalbpartropaurglab+len(mySp)
    
    # 4 spp
    if 'alb' in mySp and 'trop' in mySp and 'par' in mySp and 'aur' in mySp and not 'glab' in mySp:
        albpartropaur +=1
        genesalbpartropaur = genesalbpartropaur+len(mySp)
    if 'alb' in mySp and 'trop' in mySp and 'par' in mySp and not 'aur' in mySp and 'glab' in mySp:
        albpartropglab +=1
        genesalbpartropglab = genesalbpartropglab+len(mySp)
    if 'alb' in mySp and not 'trop' in mySp and 'par' in mySp and 'aur' in mySp and 'glab' in mySp:
        albparaurglab +=1
        genesalbparaurglab = genesalbparaurglab+len(mySp)
    if 'alb' in mySp and  'trop' in mySp and not 'par' in mySp and 'aur' in mySp and 'glab' in mySp:
        albtropaurglab +=1
        genesalbtropaurglab = genesalbtropaurglab+len(mySp)
    if not 'alb' in mySp and  'trop' in mySp and  'par' in mySp and 'aur' in mySp and 'glab' in mySp:
        partropaurglab +=1
        genespartropaurglab = genespartropaurglab+len(mySp)
    
    # 3 species
    
    if  'alb' in mySp and  'trop' in mySp and  'par' in mySp and not 'aur' in mySp and not 'glab' in mySp:
        albpartrop +=1
        genesalbpartrop = genesalbpartrop+len(mySp)
    if  'alb' in mySp and not  'trop' in mySp and  'par' in mySp and 'aur' in mySp and not 'glab' in mySp:
        albparaur +=1
        genesalbparaur = genesalbparaur+len(mySp)
    if  'alb' in mySp and not  'trop' in mySp and  'par' in mySp and not 'aur' in mySp and 'glab' in mySp:
        albparglab +=1
        genesalbparglab = genesalbparglab+len(mySp)
    if  'alb' in mySp and  'trop' in mySp and not 'par' in mySp and 'aur' in mySp and not 'glab' in mySp:
        albtropaur +=1
        genesalbtropaur = genesalbtropaur+len(mySp)
    if  'alb' in mySp and  'trop' in mySp and not 'par' in mySp and not 'aur' in mySp and 'glab' in mySp:
        albtropglab +=1
        genesalbtropglab = genesalbtropglab+len(mySp)
    if  'alb' in mySp and not 'trop' in mySp and not 'par' in mySp and  'aur' in mySp and 'glab' in mySp:
        albaurglab +=1
        genesalbaurglab = genesalbaurglab+len(mySp)
    if  not 'alb' in mySp and 'trop' in mySp and 'par' in mySp and  'aur' in mySp and not 'glab' in mySp:
        partropaur +=1
        genespartropaur = genespartropaur+len(mySp)
    if  not 'alb' in mySp and 'trop' in mySp and 'par' in mySp and not 'aur' in mySp and  'glab' in mySp:
        partropglab +=1
        genespartropglab = genespartropglab+len(mySp)
    if  not 'alb' in mySp and not 'trop' in mySp and 'par' in mySp and  'aur' in mySp and  'glab' in mySp:
        paraurglab +=1
        genesparaurglab = genesparaurglab+len(mySp)
    if  not 'alb' in mySp and 'trop' in mySp and not 'par' in mySp and  'aur' in mySp and  'glab' in mySp:
        tropaurglab +=1
        genestropaurglab = genestropaurglab+len(mySp)
        
    # 2 spp
    
    if  'alb' in mySp and not 'trop' in mySp and  'par' in mySp and not 'aur' in mySp and not 'glab' in mySp:
        albpar +=1
        genesalbpar = genesalbpar+len(mySp)
    if  'alb' in mySp and  'trop' in mySp and not 'par' in mySp and not 'aur' in mySp and not 'glab' in mySp:
        albtrop +=1
        genesalbtrop = genesalbtrop+len(mySp)
    if  'alb' in mySp and not 'trop' in mySp and not 'par' in mySp and 'aur' in mySp and not 'glab' in mySp:
        albaur +=1
        genesalbaur = genesalbaur+len(mySp)
    if  'alb' in mySp and not 'trop' in mySp and not 'par' in mySp and not 'aur' in mySp and 'glab' in mySp:
        albglab +=1
        genesalbglab = genesalbglab+len(mySp)
    if  not 'alb' in mySp and  'trop' in mySp and  'par' in mySp and not 'aur' in mySp and not 'glab' in mySp:
        partrop +=1
        genespartrop = genespartrop+len(mySp)
    if  not 'alb' in mySp and not 'trop' in mySp and  'par' in mySp and  'aur' in mySp and not 'glab' in mySp:
        paraur +=1
        genesparaur = genesparaur+len(mySp)
    if  not 'alb' in mySp and not 'trop' in mySp and  'par' in mySp and not 'aur' in mySp and  'glab' in mySp:
        parglab +=1
        genesparglab = genesparglab+len(mySp)
    if  not 'alb' in mySp and  'trop' in mySp and not 'par' in mySp and  'aur' in mySp and not 'glab' in mySp:
        tropaur +=1
        genestropaur = genestropaur+len(mySp)
    if  not 'alb' in mySp and  'trop' in mySp and not 'par' in mySp and not 'aur' in mySp and  'glab' in mySp:
        tropglab +=1
        genestropglab = genestropglab+len(mySp)
    if  not 'alb' in mySp and not 'trop' in mySp and not 'par' in mySp and 'aur' in mySp and  'glab' in mySp:
        aurglab +=1
        genesaurglab = genesaurglab+len(mySp)
        
    # 1 spp
 
    if 'alb' in mySp:
        albF +=1; albG += len(mySp)
    if 'trop' in mySp:
        tropF +=1; tropG += len(mySp)
    if 'par' in mySp:
        parF +=1; parG += len(mySp)
    if 'glab' in mySp:
        glabF +=1; glabG += len(mySp)
    if 'aur' in mySp:
		aurF +=1; aurG += len(mySp)

"""
print "albpar", albpar
print "albtrop", albtrop
print "albaur", albaur
print "albglab", albglab
print "partrop", partrop
print "paraur", paraur
print "parglab", parglab
print "tropaur", tropaur
print "tropglab", tropglab
print "aurglab", aurglab
print "albpartrop", albpartrop
print "albparaur", albparaur
print "albparglab", albparglab
print "albtropaur", albtropaur
print "albtropglab", albtropglab
print "albaurglab", albaurglab
print "partropaur", partropaur
print "partropglab", partropglab
print "paraurglab", paraurglab
print "tropaurglab", tropaurglab
print "albpartropaur", albpartropaur
print "albpartropglab", albpartropglab
print "albparaurglab", albparaurglab
print "albtropaurglab", albtropaurglab
print "partropaurglab", partropaurglab
print "albpartropaurglab", albpartropaurglab
"""


out3.write('rm(list=ls())\n')
out3.write('library(VennDiagram)\n')
out3.write('setwd("%s")\n' % os.path.abspath(os.path.dirname(sys.argv[8])))
out3.write('alb <-%s\n' % (int(calb)-(calbm -albF)))
out3.write('par <-%s\n' % (int(cpar)-(cparm-parF)))
out3.write('trop <-%s\n' % (int(ctrop)-(ctropm-tropF)))
out3.write('aur <-%s\n' % (int(caur)-(caurm-aurF)))
out3.write('glab <-%s\n' % (int(cglab)-(cglabm-glabF)))


out3.write('overlapalbpar <- %s\n' % (albpar+albpartrop+albparaur+albparglab+albpartropaur+albpartropglab+albparaurglab+albpartropaurglab))
out3.write('overlapalbtrop <- %s\n' % (albtrop+albpartrop+albtropaur+albtropglab+albpartropaur+albpartropglab+albtropaurglab+albpartropaurglab))
out3.write('overlapalbaur <- %s\n' % (albaur+albparaur+albtropaur+albaurglab+albpartropaur+albparaurglab+albtropaurglab+albpartropaurglab))
out3.write('overlapalbglab <- %s\n' % (albglab+albparglab+albtropglab+albaurglab+albpartropglab+albparaurglab+albtropaurglab+albpartropaurglab))
out3.write('overlappartrop <- %s\n' % (partrop+albpartrop+partropaur+partropglab+albpartropaur+albpartropglab+partropaurglab+albpartropaurglab))
out3.write('overlapparaur <- %s\n' % (paraur+albparaur+partropaur+paraurglab+albpartropaur+albparaurglab+partropaurglab+albpartropaurglab))
out3.write('overlapparglab <- %s\n' % (parglab+albparglab+partropglab+paraurglab+albpartropglab+albparaurglab+partropaurglab+albpartropaurglab))
out3.write('overlaptropaur <- %s\n' % (tropaur+albtropaur+partropaur+tropaurglab+albpartropaur+albtropaurglab+partropaurglab+albpartropaurglab))
out3.write('overlaptropglab <- %s\n' % (tropglab+albtropglab+partropglab+tropaurglab+albpartropglab+albtropaurglab+partropaurglab+albpartropaurglab))
out3.write('overlapaurglab <- %s\n' % (aurglab+albaurglab+paraurglab+tropaurglab+albparaurglab+albtropaurglab+partropaurglab+albpartropaurglab))
out3.write('overlapalbpartrop <- %s\n' % (albpartrop+albpartropaur+albpartropglab+albpartropaurglab))
out3.write('overlapalbparaur <- %s\n' % (albparaur+albpartropaur+albparaurglab+albpartropaurglab))
out3.write('overlapalbparglab <- %s\n' % (albparglab+albpartropglab+albparaurglab+albpartropaurglab))
out3.write('overlapalbtropaur <- %s\n' % (albtropaur+albpartropaur+albtropaurglab+albpartropaurglab))
out3.write('overlapalbtropglab <- %s\n' % (albtropglab+albpartropglab+albtropaurglab+albpartropaurglab))
out3.write('overlapalbaurglab <- %s\n' % (albaurglab+albparaurglab+albtropaurglab+albpartropaurglab))
out3.write('overlappartropaur <- %s\n' % (partropaur+albpartropaur+partropaurglab+albpartropaurglab))
out3.write('overlappartropglab <- %s\n' % (partropglab+albpartropglab+partropaurglab+albpartropaurglab))
out3.write('overlapparaurglab <- %s\n' % (paraurglab+albparaurglab+partropaurglab+albpartropaurglab))
out3.write('overlaptropaurglab <- %s\n' % (tropaurglab+albtropaurglab+partropaurglab+albpartropaurglab))
out3.write('overlapalbpartropaur <- %s\n' % (albpartropaur+albpartropaurglab))
out3.write('overlapalbpartropglab <- %s\n' % (albpartropglab+albpartropaurglab))
out3.write('overlapalbparaurglab <- %s\n' % (albparaurglab+albpartropaurglab))
out3.write('overlapalbtropaurglab <- %s\n' % (albtropaurglab+albpartropaurglab))
out3.write('overlappartropaurglab <- %s\n' % (partropaurglab+albpartropaurglab))
out3.write('overlapalbpartropaurglab <- %s\n' % (albpartropaurglab))

out3.write('venn_plot<-draw.quintuple.venn(area1=alb, area2=par, area3=trop, area4=aur,area5=glab,category = c("CALB", "CPAR", "CTROP","CAUR" ,"CGLAB"),\
cex=2, n12=overlapalbpar, n13=overlapalbtrop, n14=overlapalbaur, n15=overlapalbglab, \
n23=overlappartrop, n24=overlapparaur, n25=overlapparglab, n34=overlaptropaur, n35=overlaptropglab, n45=overlapaurglab, \
n123=overlapalbpartrop, n124=overlapalbparaur, n125=overlapalbparglab, n134=overlapalbtropaur, n135=overlapalbtropglab, n145=overlapalbaurglab, \
n234=overlappartropaur, n235=overlappartropglab, n245=overlapparaurglab,n345=overlaptropaurglab, \
n1234=overlapalbpartropaur, n1235=overlapalbpartropglab, n1245=overlapalbparaurglab, n1345=overlapalbtropaurglab, n2345=overlappartropaurglab, n12345=overlapalbpartropaurglab,\n fill = c("red", "blue", "green", "orange","pink"))\n')

out3.write('png(filename = "families.png",units="in", width=7, heigh=7, res=1200)\n')
out3.write('grid.draw(venn_plot)\n')
out3.write('dev.off()\n')
'''

out4.write('rm(list=ls())\n')
out4.write('library(VennDiagram)\n')
out4.write('setwd("%s")\n' % os.path.abspath(os.path.dirname(sys.argv[8])))
out4.write('alb <-%s\n' % ((int(calb)-calbm)+albG))
out4.write('par <-%s\n' % ((int(cpar)-cparm)+parG))
out4.write('trop <-%s\n' % ((int(ctrop)-ctropm)+tropG))
out4.write('glab <-%s\n' % ((int(cglab)-cglabm)+glabG))
out4.write('overlapalbpar <- %s\n' % (genesalbpar + genesalbpartrop + genesalbparglab + genesalbpartropglab))
out4.write('overlapalbtrop <- %s\n' % (genesalbtrop + genesalbtropglab + genesalbpartrop + genesalbpartropglab))
out4.write('overlapalbglab <- %s\n' % (genesalbglab + genesalbparglab + genesalbtropglab + genesalbpartropglab))
out4.write('overlappartrop <- %s\n' % (genespartrop + genesalbpartrop + genespartropglab + genesalbpartropglab))
out4.write('overlapparglab <- %s\n' % (genesparglab + genesalbparglab + genespartropglab + genesalbpartropglab))
out4.write('overlaptropglab <- %s\n' % (genestropglab + genesalbtropglab + genespartropglab + genesalbpartropglab))
out4.write('overlapalbpartrop <- %s\n' % (genesalbpartrop + genesalbpartropglab))
out4.write('overlapalbparglab <- %s\n' % (genesalbparglab + genesalbpartropglab))
out4.write('overlapalbtropglab <- %s\n' % (genesalbtropglab + genesalbpartropglab))
out4.write('overlappartropglab <- %s\n' % (genespartropglab + genesalbpartropglab))
out4.write('overlapalbpartropglab <- %s\n' % (genesalbpartropglab))
out4.write('venn_plot<-draw.quad.venn(area1=alb, area2=par, area3=trop, area4=glab, category = c("CALB", "CPAR", "CTROP", "CGLAB"), cex=2,\n n12=overlapalbpar, n13=overlapalbtrop, n14=overlapalbglab,\n n23=overlappartrop, n24=overlapparglab, n34=overlaptropglab,\n n123=overlapalbpartrop, n124=overlapalbparglab, n134=overlapalbtropglab,\n n234=overlappartropglab, n1234=overlapalbpartropglab,\n fill = c("red", "blue", "green", "orange"))\n')

out4.write('png(filename = "genes.png",units="in", width=7, heigh=7, res=1200)\n')
out4.write('grid.draw(venn_plot)\n')
out4.write('dev.off()\n')

'''


out4.write('rm(list=ls())\n')
out4.write('library(VennDiagram)\n')
out4.write('setwd("%s")\n' % os.path.abspath(os.path.dirname(sys.argv[8])))
out4.write('alb <-%s\n' % (int(calb)-(calbm -albG)))
out4.write('par <-%s\n' % (int(cpar)-(cparm-parG)))
out4.write('trop <-%s\n' % (int(ctrop)-(ctropm-tropG)))
out4.write('aur <-%s\n' % (int(caur)-(caurm-aurG)))
out4.write('glab <-%s\n' % (int(cglab)-(cglabm-glabG)))


out4.write('overlapalbpar <- %s\n' % (genesalbpar+genesalbpartrop+genesalbparaur+genesalbparglab+genesalbpartropaur+genesalbpartropglab+genesalbparaurglab+genesalbpartropaurglab))
out4.write('overlapalbtrop <- %s\n' % (genesalbtrop+genesalbpartrop+genesalbtropaur+genesalbtropglab+genesalbpartropaur+genesalbpartropglab+genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlapalbaur <- %s\n' % (genesalbaur+genesalbparaur+genesalbtropaur+genesalbaurglab+genesalbpartropaur+genesalbparaurglab+genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlapalbglab <- %s\n' % (genesalbglab+genesalbparglab+genesalbtropglab+genesalbaurglab+genesalbpartropglab+genesalbparaurglab+genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlappartrop <- %s\n' % (genespartrop+genesalbpartrop+genespartropaur+genespartropglab+genesalbpartropaur+genesalbpartropglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapparaur <- %s\n' % (genesparaur+genesalbparaur+genespartropaur+genesparaurglab+genesalbpartropaur+genesalbparaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapparglab <- %s\n' % (genesparglab+genesalbparglab+genespartropglab+genesparaurglab+genesalbpartropglab+genesalbparaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlaptropaur <- %s\n' % (genestropaur+genesalbtropaur+genespartropaur+genestropaurglab+genesalbpartropaur+genesalbtropaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlaptropglab <- %s\n' % (genestropglab+genesalbtropglab+genespartropglab+genestropaurglab+genesalbpartropglab+genesalbtropaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapaurglab <- %s\n' % (genesaurglab+genesalbaurglab+genesparaurglab+genestropaurglab+genesalbparaurglab+genesalbtropaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapalbpartrop <- %s\n' % (genesalbpartrop+genesalbpartropaur+genesalbpartropglab+genesalbpartropaurglab))
out4.write('overlapalbparaur <- %s\n' % (genesalbparaur+genesalbpartropaur+genesalbparaurglab+genesalbpartropaurglab))
out4.write('overlapalbparglab <- %s\n' % (genesalbparglab+genesalbpartropglab+genesalbparaurglab+genesalbpartropaurglab))
out4.write('overlapalbtropaur <- %s\n' % (genesalbtropaur+genesalbpartropaur+genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlapalbtropglab <- %s\n' % (genesalbtropglab+genesalbpartropglab+genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlapalbaurglab <- %s\n' % (genesalbaurglab+genesalbparaurglab+genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlappartropaur <- %s\n' % (genespartropaur+genesalbpartropaur+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlappartropglab <- %s\n' % (genespartropglab+genesalbpartropglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapparaurglab <- %s\n' % (genesparaurglab+genesalbparaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlaptropaurglab <- %s\n' % (genestropaurglab+genesalbtropaurglab+genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapalbpartropaur <- %s\n' % (genesalbpartropaur+genesalbpartropaurglab))
out4.write('overlapalbpartropglab <- %s\n' % (genesalbpartropglab+genesalbpartropaurglab))
out4.write('overlapalbparaurglab <- %s\n' % (genesalbparaurglab+genesalbpartropaurglab))
out4.write('overlapalbtropaurglab <- %s\n' % (genesalbtropaurglab+genesalbpartropaurglab))
out4.write('overlappartropaurglab <- %s\n' % (genespartropaurglab+genesalbpartropaurglab))
out4.write('overlapalbpartropaurglab <- %s\n' % (genesalbpartropaurglab))

out4.write('venn_plot<-draw.quintuple.venn(area1=alb, area2=par, area3=trop, area4=aur,area5=glab,category = c("CALB", "CPAR", "CTROP","CAUR" ,"CGLAB"),\
cex=2, n12=overlapalbpar, n13=overlapalbtrop, n14=overlapalbaur, n15=overlapalbglab, \
n23=overlappartrop, n24=overlapparaur, n25=overlapparglab, n34=overlaptropaur, n35=overlaptropglab, n45=overlapaurglab, \
n123=overlapalbpartrop, n124=overlapalbparaur, n125=overlapalbparglab, n134=overlapalbtropaur, n135=overlapalbtropglab, n145=overlapalbaurglab, \
n234=overlappartropaur, n235=overlappartropglab, n245=overlapparaurglab,n345=overlaptropaurglab, \
n1234=overlapalbpartropaur, n1235=overlapalbpartropglab, n1245=overlapalbparaurglab, n1345=overlapalbtropaurglab, n2345=overlappartropaurglab, n12345=overlapalbpartropaurglab,\n fill = c("red", "blue", "green", "orange","pink"))\n')

out4.write('png(filename = "genes.png",units="in", width=7, heigh=7, res=1200)\n')
out4.write('grid.draw(venn_plot)\n')
out4.write('dev.off()\n')
