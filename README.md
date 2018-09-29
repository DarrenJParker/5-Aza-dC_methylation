# 5-Aza-dC_methylation

This is the repository for the collected scripts used in the study *"5-aza-2’-deoxycytidine alters methylation patterns in insects on a genome-wide scale."* currently under preparation.

## DATA

* The locus (CpG) level file **all_methylKit_fix.csv**, contains the number of methylated and unmethylated reads per locus, can be obtained from `https://www.dropbox.com/s/s3f0amznm4tptzt/all_methylKit_fix.csv?dl=0` (Dropbox link for now, will put into long term storage later). Once downloaded please place this file into the `Data` directory. 

* Gene annotation file **all_methylKit_fix_gene_closestgene_exon_CDS_5UTR_3UTR_annot.txt**, contains gene information for each locus in the experiment. 

* Sample ID to sample infomation table **BSseq_sample_info.csv**

## Locus-level analysis

### PCA and t-SNE

* First calculate the proportion of reads that are methylated for each locus and each sample:

`python3 all_methylKit_fix_tidier.py -i ./Data/all_methylKit_fix.csv -o all_methylKit_fix`

* Then run `BS_PCA_MDS_tsne.R`

## Gene-level analysis

* First sum locus counts to gene-level counts:

`python3 all_methylKit_fix_sumlocibygene.py -a ./Data/all_methylKit_fix_gene_closestgene_exon_CDS_5UTR_3UTR_annot.txt -i ./Data/all_methylKit_fix.csv -o all_methylKit_fix -d 1000`

### PCA and t-SNE

* First calculate the proportion of reads that are methylated for each gene and each sample:

`python3 all_methylKit_fix_tidier.py -i all_methylKit_fix_gene_level.csv -o all_methylKit_fix_gene_level`

* Then run `BS_PCA_MDS_tsne_gene_level.R`

### GLM


### Additional scripts


## Infomation on running scripts

* All scripts should be run from the directory they are in. Output directories will be created to store output as the code is run. 
* All python scripts were made using python 3.5. All contain help information which can be displayed by specifying `python [name of script] -h`

## Abbreviations


