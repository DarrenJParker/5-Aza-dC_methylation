library(ggplot2)
library(cowplot)
library(tidyr)
library(cowplot)

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
# [13] "[1] stats     graphics  grDevices utils     datasets  methods   base     "                                                                                                                      
# [14] ""                                                                                                                                                                                               
# [15] "other attached packages:"                                                                                                                                                                       
# [16] "[1] tidyr_0.7.2   cowplot_0.8.0 ggplot2_2.2.1"                                                                                                                                                  
# [17] ""                                                                                                                                                                                               
# [18] "loaded via a namespace (and not attached):"                                                                                                                                                     
# [19] " [1] Rcpp_0.12.13     assertthat_0.2.0 dplyr_0.7.4      R6_2.2.2         grid_3.4.1       plyr_1.8.4       gtable_0.2.0     magrittr_1.5     scales_0.5.0     rlang_0.1.6      lazyeval_0.2.1  "
# [20] "[12] bindrcpp_0.2     glue_1.2.0       purrr_0.2.4      munsell_0.4.3    compiler_3.4.1   pkgconfig_2.0.1  colorspace_1.3-2 bindr_0.1        tidyselect_0.2.2 tibble_1.3.4    "                 


#data

dat1 <- read.csv("all_methylKit_fix_gene_level_glm_tidied_with_locID.csv", header = T)

####### filter any with NAs (this will remove genes that had 0 coverage for all samples in a treatment group)

dat2 <- dat1 %>% drop_na(
treat5aza_dC_24h,treat5aza_dC_48h,
collection24h,collection48h,
treat5aza_dC_24h_collection24h,treat5aza_dC_48h_collection24h,treat5aza_dC_24h_collection48h,treat5aza_dC_48h_collection48h,                     
treat5aza_dC_24h_tiss_typeHeads,treat5aza_dC_48h_tiss_typeHeads,
collection24h_tiss_typeHeads,collection48h_tiss_typeHeads,
treat5aza_dC_24h_collection24h_tiss_typeHeads, treat5aza_dC_48h_collection24h_tiss_typeHeads, treat5aza_dC_24h_collection48h_tiss_typeHeads, treat5aza_dC_48h_collection48h_tiss_typeHeads
)

length(dat1[,1]) - length(dat2[,1])

## drops 3183 genes

##### FDR correct

dat2$FDR_treat                            <- p.adjust(dat2$p_treat, method = "BH")
dat2$FDR_collection                       <- p.adjust(dat2$p_collection, method = "BH")
dat2$FDR_tiss_type                        <- p.adjust(dat2$p_tiss_type, method = "BH")
dat2$FDR_treat_by_collection              <- p.adjust(dat2$p_treat_by_collection, method = "BH")
dat2$FDR_treat_by_tiss_type               <- p.adjust(dat2$p_treat_by_tiss_type, method = "BH")
dat2$FDR_collection_by_tiss_type          <- p.adjust(dat2$p_collection_by_tiss_type, method = "BH")
dat2$FDR_treat_by_collection_by_tiss_type <- p.adjust(dat2$p_treat_by_collection_by_tiss_type, method = "BH")


#### N sig 

length(subset(dat2, dat2$FDR_treat < 0.05)[,1])
length(subset(dat2, dat2$FDR_collection < 0.05)[,1])
length(subset(dat2, dat2$FDR_tiss_type  < 0.05)[,1])
length(subset(dat2, dat2$FDR_treat_by_collection < 0.05)[,1])
length(subset(dat2, dat2$FDR_treat_by_tiss_type  < 0.05)[,1])
length(subset(dat2, dat2$FDR_collection_by_tiss_type  < 0.05)[,1])
length(subset(dat2, dat2$FDR_treat_by_collection_by_tiss_type < 0.05)[,1])


######### output datafile

write.csv(dat2, "all_methylKit_fix_gene_level_glm_tidied_with_locID_FDR_corrected.csv", row.names=FALSE)


######### let's plot a few things

dat2_sig_treat      <- subset(dat2, dat2$FDR_treat < 0.05)
dat2_sig_collection <- subset(dat2, dat2$FDR_collection < 0.05)
dat2_sig_tiss_type  <- subset(dat2, dat2$FDR_tiss_type < 0.05)


### treatment

treat_df_sig_treat <- as.data.frame(cbind(
c(dat2_sig_treat$treat5aza_dC_24h,dat2_sig_treat$treat5aza_dC_48h),
c(rep("treat5aza_dC_24h", length(dat2_sig_treat$treat5aza_dC_24h)), rep("treat5aza_dC_48h", length(dat2_sig_treat$treat5aza_dC_48h)))
))

colnames(treat_df_sig_treat) <- c("coeff", "treat")
treat_df_sig_treat$coeff <-  as.numeric(as.character(treat_df_sig_treat$coeff))
str(treat_df_sig_treat)



box_treat_sig_2 <- ggplot(treat_df_sig_treat, aes(treat, coeff)) + 
	theme_bw() +
	geom_boxplot(aes(fill = factor(treat)), outlier.size = 0, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	scale_fill_manual(values=c("#999999", "#999999")) + 
	ggtitle(paste("Treatment\nN =" , length(dat2_sig_treat[,1]))) +
	ylab("coeff (log odds of methylation)")
box_treat_sig_2 <- box_treat_sig_2 + theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="none")

### collection

collection_df_sig_collection <- as.data.frame(cbind(
c(dat2_sig_collection$collection24h,dat2_sig_collection$collection48h),
c(rep("collection24h", length(dat2_sig_collection$collection24h)), rep("collection48h", length(dat2_sig_collection$collection48h)))
))

colnames(collection_df_sig_collection) <- c("coeff", "collection")
collection_df_sig_collection$coeff <-  as.numeric(as.character(collection_df_sig_collection$coeff))
str(collection_df_sig_collection)


box_collection_sig_2 <- ggplot(collection_df_sig_collection, aes(collection, coeff)) + 
	theme_bw() +
	geom_boxplot(aes(fill = factor(collection)), outlier.size = 0, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	scale_fill_manual(values=c("#999999", "#999999")) + 
	ggtitle(paste("Collection\nN =" , length(dat2_sig_collection[,1]))) +
	ylab("coeff (log odds of methylation)")
box_collection_sig_2 <- box_collection_sig_2 + theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")

### tiss-type

tiss_type_df_sig_tiss_type <- as.data.frame(cbind(
c(dat2_sig_tiss_type$tiss_typeHeads),
c(rep("tiss_typeHeads", length(dat2_sig_tiss_type$tiss_typeHeads)))
))

colnames(tiss_type_df_sig_tiss_type) <- c("coeff", "tiss_type")
tiss_type_df_sig_tiss_type$coeff <-  as.numeric(as.character(tiss_type_df_sig_tiss_type$coeff))
str(tiss_type_df_sig_tiss_type)



box_tiss_type_sig_2 <- ggplot(tiss_type_df_sig_tiss_type, aes(tiss_type, coeff)) + 
	theme_bw() +
	geom_boxplot(aes(fill = factor(tiss_type)), outlier.size = 0, outlier.shape = NA) +
	coord_cartesian(ylim=c(-3,3)) +
	scale_fill_manual(values=c("#999999")) + 
	ggtitle(paste("Tissue\nN =" , length(dat2_sig_tiss_type[,1]))) +
	ylab("coeff (log odds of methylation)")
box_tiss_type_sig_2 <- box_tiss_type_sig_2 + theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.position="none")


######### output

pdf("log_odds_of_meth.pdf", width = 10, height = 5)
plot_grid(box_treat_sig_2, box_collection_sig_2, box_tiss_type_sig_2 , labels=c("A", "B", "C"), ncol = 3, nrow = 1, align = "hv")
dev.off()
getwd() ## where has my plot gone....?


########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "BSseq_analysis_of_glm_gene_level.R_sessionInfo.txt")

