This GitHub page and directory contains all scripts and data files necessary to fully reproduce the results of the paper by Hovhannisyan and Gabaldon "The lncRNA landscape of Candida yeast pathogens".

The directory contains several sub-directories for different types of analyses, which we describe below. Each sub-directory contains annotated scripts, data files and README.txt files describing the entire analysis. If analyses are done in R, the sub-directories also contain sessionInfo.txt files listing all R packages and versions.
Note that most of the analyses are computationally very demanding and, if done from start to end, the project might require ~10 terabytes of space. Hence for some of the steps we highly recommend using HPC resources. 


Below you can find general descriptions for each sub-directory in the order by which they need to be executed.

1. data_fetching_and_QC_and_strand_detection - this is the initial part of the entire workflow containing scripts and data files necessary to obtain, perform QC and trimming of all public RNA-Seq datasets used in our study. Note that downloading and processing the data is very time and space consuming since it encompasses >2600 RNA-Seq samples.
2. mapping_assembly_lncRNA_prediction - the next step of the workflow performs (among other analyses) mapping of RNA-Seq data to the corresponding reference genomes (which are located in reference_genomes_and_annotations folder), genome-guided transcriptome assembly, transcriptomes merging and coding potential assessments, which altogether constitute the process of lncRNA predictions. Again this step is computationally very expensive and ideally requires HPC resources. 
3. DE_analysis - this is the third step of the workflow (although the results of these analyses are described at the end of the paper) which performs differential expression analysis of lncRNAs across the time course of epithelial cell infection. These analyses are performed in R and can be done using an ordinary desktop computer.
4. family_classification_and_motif_analysis - this folder contains scripts and data for performing lncRNA family classifications based on BLAST reciprocal hits, secondary structure and synteny. Additionally it performs motif enrichment analysis. Structural analysis is computationally very demanding and would require an HPC and 3-4 terabytes of space.
5. repeat_calling - this folder contains scripts for perfromting repeat calling and finding the numer of lncRNAs overlapping repets. The final results of this analysis are also available in the next directory to perform final visualization. 
6. ploting_and_network_analysis - the final step of the workflow is used to perform most of the plotting tasks of the entire project, along with functional and network analysis. This step is done in R and can be executed in a normal desktop (although the network analysis might be quite slow).

Each of the above described sub-directories has its own README.txt file with more detailed descriptions on how to execute the scripts and which software are necessary.
Note that the entire workflow uses a lot of different software (and obviously many more dependencies), so we would recommend to install and handle them using package managers such as conda.

Also note that for the data analysis convenience we use shortened names of the studied species throughout the workflow, which appear in file/folder names. Those are:
calb - C. albicans;
ctrop - C. tropicalis;
cpar - C. parapsilosis;
caur - C. auris;
cglab - C. glabrata; 

If there are any questions or requests, please email them to grant.hovhannisyan@gmail.com.   
