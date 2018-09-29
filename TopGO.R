### topGO
# install
# source("http://bioconductor.org/biocLite.R") 
# biocLite() 
# source("http://bioconductor.org/biocLite.R")   
# biocLite("topGO")
# biocLite("ALL")
# biocLite("affyLib")

library(topGO)
library(ALL)
capture.output(sessionInfo())

# [1] "R version 3.4.1 (2017-06-30)"                                                                                                                                                                      
 # [2] "Platform: x86_64-apple-darwin15.6.0 (64-bit)"                                                                                                                                                      
 # [3] "Running under: macOS High Sierra 10.13.6"                                                                                                                                                          
 # [4] ""                                                                                                                                                                                                  
 # [5] "Matrix products: default"                                                                                                                                                                          
 # [6] "BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib"                                                                                                                 
 # [7] "LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib"                                                                                                               
 # [8] ""                                                                                                                                                                                                  
 # [9] "locale:"                                                                                                                                                                                           
# [10] "[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8"                                                                                                                                 
# [11] ""                                                                                                                                                                                                  
# [12] "attached base packages:"                                                                                                                                                                           
# [13] "[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     "                                                                                                     
# [14] ""                                                                                                                                                                                                  
# [15] "other attached packages:"                                                                                                                                                                          
# [16] " [1] ALL_1.18.0           topGO_2.28.0         SparseM_1.77         GO.db_3.4.1          AnnotationDbi_1.38.2 IRanges_2.10.5       S4Vectors_0.14.7     Biobase_2.36.2       graph_1.54.0        " 
# [17] "[10] BiocGenerics_0.22.1 "                                                                                                                                                                         
# [18] ""                                                                                                                                                                                                  
# [19] "loaded via a namespace (and not attached):"                                                                                                                                                        
# [20] " [1] Rcpp_0.12.13       matrixStats_0.52.2 lattice_0.20-35    digest_0.6.12      grid_3.4.1         DBI_0.7            RSQLite_2.0        rlang_0.1.6        blob_1.1.0         bit64_0.9-7       "
# [21] "[11] bit_1.1-12         compiler_3.4.1     pkgconfig_2.0.1    memoise_1.1.0      tibble_1.3.4      "                                                                                               


#### load annotations

geneID2GO_Nvit           <- readMappings(file = "Data/Nvit_1.2_GOs_from_hymenopteramine_270918.tsv_fortopgo.txt") # 
geneID2GO_Dmel_union     <- readMappings(file = "Data/Drosp_GOs_from_hymenopteramine_for_Nvit_1.2_orths_270918.tsv_fortopgo_union.txt") #  D. melanogaster GO-terms | GOs from multiple orths combined together
geneID2GO_Dmel_intersect <- readMappings(file = "Data/Drosp_GOs_from_hymenopteramine_for_Nvit_1.2_orths_270918.tsv_fortopgo_intersection.txt") #  D. melanogaster GO-terms | only GO-terms that were shared between all orthologs were kept

###############################################################################################################################################
#### read in tables with genename and sig value

all_meth_data <- read.csv("all_methylKit_fix_gene_level_glm_tidied_with_locID_FDR_corrected.csv")

#### need to get the genes as a named numeric vector

treat_list        <- all_meth_data$FDR_treat
names(treat_list) <- all_meth_data$gene_name

head(treat_list, n = 10)

## set the sig_for_go to set the threshold for significant GO terms from the output of the GSEA
run_enrichment <- function(genelist, ref, sig_for_GO){
	
	### make rule for classing sig / non-sig - note this rule is not used for the GSEA
	
	topDiffGenes <- function(allScore) {return(allScore < 0.05)}
	# topDiffGenes <- function(allScore) {return(allScore < 1)} ## as a check - setting to one gives the same pvalues for the GSEA
	
	#### make GOdata object
	#### setting node size as 5 so at least 5 genes must be annot per GO terms 
	#### do enrichment test
	
	GODATA_BP = new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 5)
	GODATA_MF = new("topGOdata", ontology = "MF", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 5)	
	GODATA_CC = new("topGOdata", ontology = "CC", allGenes = genelist, geneSel = topDiffGenes,  annot = annFUN.gene2GO, gene2GO = ref, nodeSize = 5)	
	
	### get N GOs used

	### get N GOs used

	GO_term_use_BP_list = GODATA_BP@graph@nodes
	GO_term_use_MF_list = GODATA_MF@graph@nodes
	GO_term_use_CC_list = GODATA_CC@graph@nodes
	GO_term_use_ALL_list = c(GO_term_use_BP_list,GO_term_use_MF_list,GO_term_use_CC_list)
	
	N_GO_term_use_BP = length(GODATA_BP@graph@nodes)
	N_GO_term_use_MF = length(GODATA_MF@graph@nodes)
	N_GO_term_use_CC = length(GODATA_CC@graph@nodes)
	N_GO_term_use_ALL = sum(N_GO_term_use_BP,N_GO_term_use_MF,N_GO_term_use_CC)
	

	result_GSEA_BP     <- runTest(GODATA_BP, algorithm = "elim", statistic = "ks")
	result_GSEA_MF     <- runTest(GODATA_MF, algorithm = "elim", statistic = "ks")
	result_GSEA_CC     <- runTest(GODATA_CC, algorithm = "elim", statistic = "ks")
	
	### combined tables
	allRes1_BP <- GenTable(GODATA_BP, GSEA = result_GSEA_BP, ranksOf = "GSEA", topNodes = length(GODATA_BP@graph@nodes), numChar = 200)
	sig_GSEA_BP_GO     = subset(allRes1_BP, allRes1_BP$GSEA < sig_for_GO)$GO.ID
	allRes1_MF <- GenTable(GODATA_MF, GSEA = result_GSEA_MF, ranksOf = "GSEA", topNodes = length(GODATA_MF@graph@nodes), numChar = 200)
	sig_GSEA_MF_GO     = subset(allRes1_MF, allRes1_MF$GSEA < sig_for_GO)$GO.ID
	allRes1_CC <- GenTable(GODATA_CC, GSEA = result_GSEA_CC, ranksOf = "GSEA", topNodes = length(GODATA_CC@graph@nodes), numChar = 200)
	sig_GSEA_CC_GO     = subset(allRes1_CC, allRes1_CC$GSEA < sig_for_GO)$GO.ID
	
	allRes1_ALL = rbind(allRes1_BP,allRes1_MF, allRes1_CC)
	sig_GSEA_ALL_GO     = c(sig_GSEA_BP_GO,sig_GSEA_MF_GO,sig_GSEA_CC_GO)
	
	## return everything!
	out_list = list("N_GO_term_use_BP" = N_GO_term_use_BP, "N_GO_term_use_MF" = N_GO_term_use_MF, "N_GO_term_use_CC" = N_GO_term_use_CC, "N_GO_term_use_ALL" = N_GO_term_use_ALL, 
	                "GO_term_use_BP_list" = GO_term_use_BP_list, "GO_term_use_MF_list" = GO_term_use_MF_list, "GO_term_use_CC_list" = GO_term_use_CC_list, "GO_term_use_ALL_list" = GO_term_use_ALL_list, 
	                "allRes1_BP" = allRes1_BP, "allRes1_MF" = allRes1_MF, "allRes1_CC" = allRes1_CC, "allRes1_ALL" = allRes1_ALL, 
	                "sig_GSEA_BP_GO" = sig_GSEA_BP_GO, "sig_GSEA_MF_GO" = sig_GSEA_MF_GO, "sig_GSEA_CC_GO" = sig_GSEA_CC_GO, "sig_GSEA_ALL_GO" = sig_GSEA_ALL_GO,
	                "GODATA_BP" = GODATA_BP, "GODATA_MF" = GODATA_MF, "GODATA_CC" = GODATA_CC) 
	return(out_list)

}

#### run the enrichment stuff (0.05)

### nr

RGL_Nvit_treat_enrich  <- run_enrichment(treat_list, geneID2GO_Nvit, 0.05)
RGL_Droso_union_treat_enrich  <- run_enrichment(treat_list, geneID2GO_Dmel_union, 0.05)
RGL_Droso_inter_treat_enrich  <- run_enrichment(treat_list, geneID2GO_Dmel_intersect, 0.05)

head(RGL_Nvit_treat_enrich$allRes1_BP, n = 100)
head(RGL_Nvit_treat_enrich$allRes1_MF, n = 100)
head(RGL_Nvit_treat_enrich$allRes1_CC, n = 100)

head(RGL_Droso_union_treat_enrich$allRes1_BP, n = 100)
head(RGL_Droso_union_treat_enrich$allRes1_MF, n = 100)
head(RGL_Droso_union_treat_enrich$allRes1_CC, n = 100)

head(RGL_Droso_inter_treat_enrich$allRes1_BP, n = 100)
head(RGL_Droso_inter_treat_enrich$allRes1_MF, n = 100)
head(RGL_Droso_inter_treat_enrich$allRes1_CC, n = 100)

### all GO out

#### Nvit
RGL_Nvit_treat_enrich_allRes1_ALL <- rbind(
RGL_Nvit_treat_enrich$allRes1_BP, 
RGL_Nvit_treat_enrich$allRes1_MF,
RGL_Nvit_treat_enrich$allRes1_CC
)
RGL_Nvit_treat_enrich_allRes1_ALL_2 <- as.data.frame(cbind(
RGL_Nvit_treat_enrich_allRes1_ALL$GO.ID, 
RGL_Nvit_treat_enrich_allRes1_ALL$Term,
c(rep("BP", length(RGL_Nvit_treat_enrich$allRes1_BP[,1])), rep("MF", length(RGL_Nvit_treat_enrich$allRes1_MF[,1])),  rep("CC", length(RGL_Nvit_treat_enrich$allRes1_CC[,1]))), 
RGL_Nvit_treat_enrich_allRes1_ALL$GSEA
))

colnames(RGL_Nvit_treat_enrich_allRes1_ALL_2) <- c("GO.ID", "Term_desc", "GO_type","GSEA_p")
head(RGL_Nvit_treat_enrich_allRes1_ALL_2)


#### Droso_union
RGL_Droso_union_treat_enrich_allRes1_ALL <- rbind(
RGL_Droso_union_treat_enrich$allRes1_BP, 
RGL_Droso_union_treat_enrich$allRes1_MF,
RGL_Droso_union_treat_enrich$allRes1_CC
)
RGL_Droso_union_treat_enrich_allRes1_ALL_2 <- as.data.frame(cbind(
RGL_Droso_union_treat_enrich_allRes1_ALL$GO.ID, 
RGL_Droso_union_treat_enrich_allRes1_ALL$Term,
c(rep("BP", length(RGL_Droso_union_treat_enrich$allRes1_BP[,1])), rep("MF", length(RGL_Droso_union_treat_enrich$allRes1_MF[,1])),  rep("CC", length(RGL_Droso_union_treat_enrich$allRes1_CC[,1]))), 
RGL_Droso_union_treat_enrich_allRes1_ALL$GSEA
))

colnames(RGL_Droso_union_treat_enrich_allRes1_ALL_2) <- c("GO.ID", "Term_desc", "GO_type","GSEA_p")
head(RGL_Droso_union_treat_enrich_allRes1_ALL_2)


#### Droso_inter
RGL_Droso_inter_treat_enrich_allRes1_ALL <- rbind(
RGL_Droso_inter_treat_enrich$allRes1_BP, 
RGL_Droso_inter_treat_enrich$allRes1_MF,
RGL_Droso_inter_treat_enrich$allRes1_CC
)
RGL_Droso_inter_treat_enrich_allRes1_ALL_2 <- as.data.frame(cbind(
RGL_Droso_inter_treat_enrich_allRes1_ALL$GO.ID, 
RGL_Droso_inter_treat_enrich_allRes1_ALL$Term,
c(rep("BP", length(RGL_Droso_inter_treat_enrich$allRes1_BP[,1])), rep("MF", length(RGL_Droso_inter_treat_enrich$allRes1_MF[,1])),  rep("CC", length(RGL_Droso_inter_treat_enrich$allRes1_CC[,1]))), 
RGL_Droso_inter_treat_enrich_allRes1_ALL$GSEA
))

colnames(RGL_Droso_inter_treat_enrich_allRes1_ALL_2) <- c("GO.ID", "Term_desc", "GO_type","GSEA_p")
head(RGL_Droso_inter_treat_enrich_allRes1_ALL_2)


##### output

write.table(RGL_Nvit_treat_enrich_allRes1_ALL_2,        "RGL_Nvit_treat_GO_terms.txt",        sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Droso_union_treat_enrich_allRes1_ALL_2, "RGL_Droso_union_treat_GO_terms.txt", sep = '\t', quote = FALSE, row.names = FALSE)
write.table(RGL_Droso_inter_treat_enrich_allRes1_ALL_2, "RGL_Droso_inter_treat_GO_terms.txt", sep = '\t', quote = FALSE, row.names = FALSE)


########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "TopGO.R_sessionInfo.txt")





