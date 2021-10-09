### This script will run predict_lncRNAs.py for all species. The last argument are files containing the stringtie assemblies of all sample per species. For example, in case of C. auris, the file would look like:

#full/path/mapping_assembly_lncRNA_prediction/caurstringtie_assembly/SRR6900282.gtf
#full/path/mapping_assembly_lncRNA_prediction/caurstringtie_assembly/SRR6900283.gtf
#full/path/mapping_assembly_lncRNA_prediction/caurstringtie_assembly/SRR6900284.gtf
#full/path/mapping_assembly_lncRNA_prediction/caurstringtie_assembly/SRR6900285.gtf
#etc

### NOTE: You need to uncomment functions in predict_lncRNAs.py to execute them. We recommend to run each or several functions step by step to follow the workflow of the script easier.

python predict_lncRNAs.py calb calb_assemblies.txt
python predict_lncRNAs.py ctrop ctrop_assemblies.txt
python predict_lncRNAs.py cpar cpar_assemblies.txt
python predict_lncRNAs.py caur caur_assemblies.txt
python predict_lncRNAs.py cglab cglab_assemblies.txt

