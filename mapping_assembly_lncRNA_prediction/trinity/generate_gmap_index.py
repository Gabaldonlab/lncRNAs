import subprocess
import os


### Before all the analyses, run the fuction below to generate gmap genome indices.
def generate_gmap_index():
    for spp in ["calb","ctrop","cpar","caur","cglab"]:
        index = "source ~/miniconda2/etc/profile.d/conda.sh && conda activate trinity && \
        gmap_build -d %s ../../reference_genomes_and_annotations/%s.fa  -D ./%s \
         && conda deactivate"%(spp,spp,spp,spp)
        print(index)
        subprocess.call(index,shell=True)

#generate_gmap_index()
