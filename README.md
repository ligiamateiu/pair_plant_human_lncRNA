# pair_plant_human_lncRNA

## Contents of the Repository:

### Input Files: 
assembly_summary.txt downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/

### Short Bash Scripts: 

`prepare_gtfs.sh`: download and extract all plant lncRNAs coordinates.

`alignment_plants10k_hg38.sh`:  download refseq FASTA files, convert them into 10-kilomer (10k-mer) fasta sequences, align these 10k-mers to the human genome (hg38), and extract the coordinates of the mapped 10k-mers that overlap with human and plant lncRNAs.

`blastn_plantlncrna_humanlncrna.sh`: extract human lncRNA sequences and create a blastn database for them, extract plant lncRNA sequences, and find sequence similarity for the paired plant - human lncRNA using blastn (specifically, megablast)

### Results:

`bams.zip`: plant 10k-mers mapping to hg38

`beds.zip`: bed files for plant 10k-mers mapping to hg38 summarized in the `tables all_plant10k_vs_human_genes.bed` and `all_plant10k_vs_human_lncRNA.bed`

`blastoutput_plant_human_lncRNA.zip`: blastn output for paired plant - human lncRNA; summarized in `table plant_lncrna_human_lncrna_1000overlap.tsv` 

![Workflow](https://github.com/ligiamateiu/pair_plant_human_lncRNA/blob/main/results/Figure1.png)
