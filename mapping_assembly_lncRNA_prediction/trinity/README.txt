This folder contains scripts to produce Trinity assemblies and compare them with Stringtie assemblies.
First, run the script generate_gmap_index.py to generate genome indices for gmap.
Next, for each species, run commands in files <species>_commands.txt. We recommend to use HPC. 
Note that these commands must be run only when all previous steps of the pipeline are finished. Also note that paths in these commands are internal for our file system, and must be changed to fit yours.

When the all the commands are succesfully finished, run compare_trinity.py. The final files for generating plots are available in this folder.
Those are:
unique_transcirpts_trinity_vs_stringtie.txt
unique_transcirpts_trinity_vs_stringtie_ux.txt
<species>/sensitivity_precision_stringtie.txt
<species>/sensitivity_precision_trinity.txt
<species>/all_class_codes_trinity.txt
<species>/all_class_codes_trinity.txt
