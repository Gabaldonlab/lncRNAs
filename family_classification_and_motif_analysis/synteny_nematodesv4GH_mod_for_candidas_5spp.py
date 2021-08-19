# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 16:52:45 2016
@author: cinta
modified by Hrant Hovhnnisyan
"""

import sys,re

def recodeList(file1, file2, string):#gene order list for spX, orthology file, ex: 'calb'
    #returns a list with the gene order of a given sp were geneIDs were renamed to myID 
    #to allow direct comparisson between species;
    #genes with no homology are also included!
    #ex: ['myID1', 'WBGene000xxx' 'myID2', 'myID3', 'myID4', 'myID5','XLOC_001892', 'myID6', 'myID7', ...]  
    mylist = []
    matched=0; nonmatched=0; found=0; lncRNA=0
    for row1 in open(file1, 'r'):  
        row1 = row1.rstrip()
        # ~ if re.search(r'MSTRG', row1):
        if row1.startswith("MSTRG"):
            #print(row1)
            mylist.append(row1)
            found=1
            lncRNA=lncRNA+1
        for row2 in open(file2 , 'r').readlines():  
            row2= row2.rstrip().split("\t")
            if string=='calb':
                if row1== row2[0]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='ctrop':
                if row1== row2[1]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='cpar':
                if row1== row2[2]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            ##    
            if string=='caur':
                if row1== row2[3]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            ###
            if string=='cglab':
                if row1== row2[4]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
                    
        if found==0:
            nonmatched=nonmatched+1
            if nonMatch=='yes':
                mylist.append(row1)
        found=0
    print "number of matched IDs for",string,": " ,matched
    print "number of NON matched IDs for",string,": " ,nonmatched
    print "number of candidate lncRNA",string,": " ,lncRNA,"\n" 
    return mylist

def dictionaryOfClusters(myidx,mylist):#myidx=positions in myList for each lncRNA; myList=renamed gene order
    #to create a dict; for each lncRNA (key) stores the number of nearby renamed geneID indicated by the user  
    mydict={}
    for idx in myidx:
        key=mylist[idx]
        val={'left':[], 'right':[], 'all':[]}
        if not key in mydict:
            mydict[key]=val
        i=1
        while i <=genesNearby:
            try:
                if idx-i >=0:                 
                    mydict[key]['left'].append(mylist[idx-i])
                    mydict[key]['all'].append(mylist[idx-i])
                mydict[key]['right'].append(mylist[idx+i])                                   
                mydict[key]['all'].append(mylist[idx+i])
            except IndexError:
                pass
            i=i+1
            
    return mydict
    
def comparingDict(sp1, sp2):
    #compares dictionaries from dictionaryOfClusters() for two species; 
    #if the number of shared genes is >= minOverlap, lncRNA are stored in myHomologs list
    if sp1=='calb':
        dict1= dictionaryOfClusters(calb_idx,calbList)
    if sp1=='ctrop':
        dict1= dictionaryOfClusters(ctrop_idx,ctropList)
    if sp1=='cpar':
        dict1= dictionaryOfClusters(cpar_idx,cparList)
    ##
    if sp1=='caur':
		dict1 = dictionaryOfClusters(caur_idx,caurList)
    ##
    if sp1=='cglab':
        dict1= dictionaryOfClusters(cglab_idx,cglabList)
        
        
    if sp2=='calb':
        dict2= dictionaryOfClusters(calb_idx,calbList)
    if sp2=='ctrop':
        dict2= dictionaryOfClusters(ctrop_idx,ctropList)
    if sp2=='cpar':
        dict2= dictionaryOfClusters(cpar_idx,cparList)
    ##    
    if sp2=='caur':
		dict2=dictionaryOfClusters(caur_idx,caurList)
		
	##
    if sp2=='cglab':
        dict2= dictionaryOfClusters(cglab_idx,cglabList)
    myHomologs=[]
    homologFound='false'
    
    
    for key1, val1 in dict1.iteritems():
        #print key1, val1
        for key2, val2 in dict2.iteritems():
            if len(set(dict1[key1]['all']).intersection(dict2[key2]['all'])) >=minOverlap:
            #at least one of the two lncRNA share genes in the left and the right side    
            #forces that both lncRNA share genes in the left and the right side    
                if len(set(dict1[key1]['right']).intersection(dict2[key2]['right'])) >=minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['left'])) >=minSideOverlap: 
                        homologFound='true'
                if len(set(dict1[key1]['right']).intersection(dict2[key2]['left'])) >=minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['right'])) >=minSideOverlap:
                        homologFound='true'
                if homologFound=='true':             
                        mytup=(key1,key2)
                        myHomologs.append(mytup)
                        homologFound='false'
    return myHomologs

#####################################

in1= sys.argv[1] #gene order list for sp1
in2= sys.argv[2] #gene order list for sp2
in3= sys.argv[3] #gene order list for sp3
in4= sys.argv[4] #gene order list for sp4
in5= sys.argv[5] #gene order list for sp5
in6= sys.argv[6] #orthology file
out= open(sys.argv[7] , 'w') #out file
genesNearby= int(sys.argv[8])
minOverlap= int(sys.argv[9])
minSideOverlap= int(sys.argv[10])
nonMatch = sys.argv[11]
temp= open('temp', 'w')

calbList=recodeList(in1, in6, 'calb')
ctropList=recodeList(in2, in6, 'ctrop')
cparList=recodeList(in3, in6, 'cpar')
caurList=recodeList(in4, in6, 'caur')
cglabList=recodeList(in5, in6, 'cglab')

for x in calbList:
    temp.write("%s\tcalb\n" %(x))       
for x in ctropList :
    temp.write("%s\tctrop\n" %(x))    
for x in cparList:
    temp.write("%s\tcpar\n" %(x))
for x in caurList:
    temp.write("%s\tcaur\n" %(x))    
for x in cglabList:
    temp.write("%s\tcglab\n" %(x))

calb_idx = [i for i, item in enumerate(calbList) if item.startswith('MSTRG')]
#list containing the positions in calbList for each lncRNA
ctrop_idx = [i for i, item in enumerate(ctropList) if item.startswith('MSTRG')]
cpar_idx = [i for i, item in enumerate(cparList) if item.startswith('MSTRG')]
caur_idx = [i for i, item in enumerate(caurList) if item.startswith('MSTRG')]
cglab_idx = [i for i, item in enumerate(cglabList) if item.startswith('MSTRG')]


for x in comparingDict('calb','ctrop'):
    out.write('calb\tctrop\t%s\n'% ('\t'.join(x)))  
for x in comparingDict('calb','cpar'):
    out.write('calb\tcpar\t%s\n'% ('\t'.join(x)))
for x in comparingDict('calb','caur'):
    out.write('calb\tcaur\t%s\n'% ('\t'.join(x)))
for x in comparingDict('calb','cglab'):
    out.write('calb\tcglab\t%s\n'% ('\t'.join(x)))
for x in comparingDict('ctrop','cpar'):
    out.write('ctrop\tcpar\t%s\n'% ('\t'.join(x)))
for x in comparingDict('ctrop','caur'):
    out.write('ctrop\tcaur\t%s\n'% ('\t'.join(x)))
for x in comparingDict('ctrop','cglab'):
    out.write('ctrop\tcglab\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cpar','caur'):     
    out.write('cpar\tcaur\t%s\n'% ('\t'.join(x)))    
for x in comparingDict('cpar','cglab'):     
    out.write('cpar\tcglab\t%s\n'% ('\t'.join(x)))
for x in comparingDict('caur','cglab'):     
    out.write('caur\tcglab\t%s\n'% ('\t'.join(x)))
    
    
for x in comparingDict('ctrop','calb'):
    out.write('ctrop\tcalb\t%s\n'% ('\t'.join(x)))  
for x in comparingDict('cpar','calb'):
    out.write('cpar\tcalb\t%s\n'% ('\t'.join(x)))
for x in comparingDict('caur','calb'):
    out.write('caur\tcalb\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cglab','calb'):
    out.write('cglab\tcalb\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cpar','ctrop'):
    out.write('cpar\tctrop\t%s\n'% ('\t'.join(x)))
for x in comparingDict('caur','ctrop'):
    out.write('caur\tctrop\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cglab','ctrop'):
    out.write('cglab\tctrop\t%s\n'% ('\t'.join(x)))
for x in comparingDict('caur','cpar'):     
    out.write('caur\tcpar\t%s\n'% ('\t'.join(x)))    
for x in comparingDict('cglab','cpar'):     
    out.write('cglab\tcpar\t%s\n'% ('\t'.join(x)))
for x in comparingDict('cglab','caur'):     
    out.write('cglab\tcaur\t%s\n'% ('\t'.join(x)))
    
