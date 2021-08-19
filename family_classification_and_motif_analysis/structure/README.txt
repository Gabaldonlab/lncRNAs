To generate secondary-stracture-based classification of lncRNAs just run

python run_rnafold_and_beagle.py

Please note that this analysis is very cumbersome and might require an HPC to run (specifically all pairwise commands of Beagle, see notes inside the script).


The folder ./reshuffled has a script reshuffle_and_run_rnafold_and_beagle.py to perform the same analysis but with random sequences. This analysis is even more 
complicated and requires a large HPC and a lot of storage. Follow the comments inside the script for details.
