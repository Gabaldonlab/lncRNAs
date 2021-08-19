#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 21:13:16 2020

@author: hrant
"""

import pandas as pd
import os
import sys

beagle_list=[]
with open(sys.argv[1], "r+") as beagle_data:
    for line in beagle_data:
        line = line.rstrip().split("|")
        id_short=line[0][1:]+"_"+line[5]
        id_1=line[0][1:] + "|" +\
            line[1] + "|" +\
            line[2] + "|" +\
            line[3] + "|" +\
            line[4] + "|" +\
            line[5]
        id_2=line[6] + "|" +\
            line[7] + "|" +\
            line[8] + "|" +\
            line[9] + "|" +\
            line[10] + "|" +\
            line[11]
            
        zscore=line[-1].split(":")[1]
        pval=line[-2].split(":")[1]

                
        #print(line[0][1:])
        beagle_list.append([id_short,id_1,id_2,float(zscore),float(pval)])

#print(beagle_dict)
        
df=pd.DataFrame(beagle_list)   
df.columns=["id_short","ids_1","id_2","zscore","pval"] 

df['zscore'] = df['zscore'].astype(float)
df['pval'] = df['pval'].astype(float)

df=(df[df["pval"]<0.01])
df=(df[df["zscore"]>3])
        
inds = df.groupby(['id_short'])['zscore'].transform(max) == df['zscore']
df = df[inds]
df.reset_index(drop=True, inplace=True)
print(df)                   

results_list = df.values.tolist()
#print (results_list)



with open("%s/results_for_clustering_script_%s.txt"%(sys.argv[3],sys.argv[2]),"w") as output:
    for el in results_list:
        #print(el[1].split("|")[-1],el[2].split("|")[-1],el[1],el[2])
        output.write("%s\t%s\t%s\t%s\n"%(el[1].split("|")[-1],el[2].split("|")[-1],el[1],el[2]))

with open("%s/results_for_clustering_script_%s.txt"%(sys.argv[3],sys.argv[2]),"a+") as output:
    for el in results_list:
    	output.write("%s\t%s\t%s\t%s\n"%(el[2].split("|")[-1],el[1].split("|")[-1],el[2],el[1]))

