To search for motifs in syntenic families, simply run

python motif_analysis.py

You can increase the number of simulations, see the comments in the script.

The script motif_analysis.py performs motif enrichment using MEME software. 

An alternative software is lncLOOM, and the analysis with that software can be done with 
cd lncLOOM
python lncLOOM.py

Again the number of simulations can be increased.

After that step is finished, run
python parse_lncLOOM.py
