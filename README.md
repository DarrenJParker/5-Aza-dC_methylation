# 5-Aza-dC_methylation

This is the repository for the collected scripts used in the study:
>Cook, N*., Parker, D. J*., Turner, F., Tauber, E., Pannebakker, B. A., Shuker, D. M. 2018. [Genome-wide disruption of DNA methylation by 5-aza-2’-deoxycytidine in the parasitoid wasp Nasonia vitripennis](https://www.biorxiv.org/content/early/2018/10/09/437202). bioRxiv.

currently under review.

## DATA

In the `Data` directory:

* The locus (CpG) level file  **GSE125388_PropMethDataAllSamples.csv.gz** **all_methylKit_fix.csv**, contains the number of methylated and unmethylated reads per locus, can be obtained from [GEO (accession number: GSE125388)](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE125388&format=file&file=GSE125388%5FPropMethDataAllSamples%2Ecsv%2Egz). Once downloaded please unzip this file and place into the `Data` directory. 

* gff file **Nasonia_vitripennis.Nvit_2.1.40.gff3**, obtained from ensembl: `ftp://ftp.ensemblgenomes.org/pub/release-40/metazoa/gff3/nasonia_vitripennis/Nasonia_vitripennis.Nvit_2.1.40.gff3.gz`. Once downloaded please unzip and place this file into the `Data` directory.

* The processed gene annotation file **all_methylKit_fix_gene_closestgene_exon_CDS_5UTR_3UTR_annot.txt**, contains gene information for each locus in the experiment. Needs to be unzipped before it can be used. 

* Sample ID to sample infomation table **BSseq_sample_info.csv**

* Refseq_ids **Nvit_OGSv1.2_official_id_map.txt** 
    * available from http://hymenopteragenome.org/nasonia/nasonia_genome_consortium/data/Nvit_OGSv1.2_official_id_map.txt

* GO-terms files from http://hymenopteragenome.org

    * **Nvit_1.2_GOs_from_hymenopteramine_270918.tsv** - GO-terms for *N. vitripennis*
    * **Drosp_GOs_from_hymenopteramine_for_Nvit_1.2_orths_270918.tsv** - GO-terms for *D. melanogaster* orthologs
    * **Gene_Orthologues_for_Nvit_1.2_from_hymenopteramine_270918.tsv** - Table for converting between *N. vitripennis* and *D. melanogaster* orthologs

### DATA - expected output

In the `Data/expected_output` directory:

* **all_methylKit_fix_gene_level_glm_tidied_with_locID_with_posinfo.csv**. Expected output file for the gene level analyses (see below), provided for convenience. 

## PCA and t-SNE

**Gene Level** 

* First sum locus counts to gene-level counts:

```
python3 all_methylKit_fix_sumlocibygene.py -a ./Data/all_methylKit_fix_gene_closestgene_exon_CDS_5UTR_3UTR_annot.txt -i ./Data/GSE125388_PropMethDataAllSamples.csv -o all_methylKit_fix -d 1000
```

* Calculate the proportion of reads that are methylated for each gene and each sample:

```
python3 all_methylKit_fix_tidier.py -i all_methylKit_fix_gene_level.csv -o all_methylKit_fix_gene_level
```

* Then run `BS_PCA_MDS_tsne_gene_level.R`

**Locus-level**

* Calculate the proportion of reads that are methylated for each locus and each sample:

```
python3 all_methylKit_fix_tidier.py -i ./Data/GSE125388_PropMethDataAllSamples.csv -o all_methylKit_fix
```

* Then run `BS_PCA_MDS_tsne.R`


## GLMs

* Firstly need to prepare file for R (one per gene):

```
python3 all_methylKit_fix_prepforglm.py -i  all_methylKit_fix_gene_level.csv -s ./Data/BSseq_sample_info.csv -o all_methylKit_fix_gene_level
```

* Then run a glm on each gene using BSseq_binomial_glm.R in a loop. Here I used a simple bash loop:

```
for i in ./single_files_for_Rglm/*.csv; do
	echo $i
	Rscript BSseq_binomial_glm.R $i
done
```

* then stick output together, add refseq_ids and position infomation.

```
python3 BSseq_binomial_glm_tidier.py -d ./single_files_for_Rglm/ -o all_methylKit_fix_gene_level
python3 add_NB_gene_info.py -t ./Data/Nvit_OGSv1.2_official_id_map.txt -i all_methylKit_fix_gene_level_glm_tidied.csv
python3 add_location_gene_info.py -g ./Data/Nasonia_vitripennis.Nvit_2.1.40_gene.gff3 -i all_methylKit_fix_gene_level_glm_tidied_with_locID.csv 
```

* note the output of these scripts should produce a file called **all_methylKit_fix_gene_level_glm_tidied_with_locID_with_posinfo.csv**. For convenience this has been added to ./Data/expected_output/

* Filter glm and correct for multiple tests with `BSseq_analysis_of_glm_gene_level.R`


## GO-term analysis

* First convert hymenoptera mine files for use in TopGO:

```
python3 hymenopteramine_GOs_to_topGO.py \
-g ./Data/Nvit_1.2_GOs_from_hymenopteramine_270918.tsv \
-l ./Data/Nasonia_vitripennis.Nvit_2.1_genenames.txt \
-d ./Data/Drosp_GOs_from_hymenopteramine_for_Nvit_1.2_orths_270918.tsv \
-o ./Data/Gene_Orthologues_for_Nvit_1.2_from_hymenopteramine_270918.tsv
```

* Then perform enrichment analyses with `TopGO.R`

## Infomation on running scripts

* All scripts should be run from the directory they are in. Output directories will be created to store output as the code is run. 
* All python scripts were made using python 3.5. All contain help information which can be displayed by specifying `python [name of script] -h`



