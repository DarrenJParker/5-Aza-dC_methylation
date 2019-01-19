library(ggplot2)
library(cowplot)
library(tidyr)
library(cowplot)
library(scales)


capture.output(sessionInfo())
 # [1] "R version 3.5.1 (2018-07-02)"                                                                                                                                                                                                     
 # [2] "Platform: x86_64-apple-darwin15.6.0 (64-bit)"                                                                                                                                                                                     
 # [3] "Running under: macOS  10.14.2"                                                                                                                                                                                                    
 # [4] ""                                                                                                                                                                                                                                 
 # [5] "Matrix products: default"                                                                                                                                                                                                         
 # [6] "BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib"                                                                                                                                                
 # [7] "LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib"                                                                                                                                              
 # [8] ""                                                                                                                                                                                                                                 
 # [9] "locale:"                                                                                                                                                                                                                          
# [10] "[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8"                                                                                                                                                                
# [11] ""                                                                                                                                                                                                                                 
# [12] "attached base packages:"                                                                                                                                                                                                          
# [13] "[1] stats     graphics  grDevices utils     datasets  methods   base     "                                                                                                                                                        
# [14] ""                                                                                                                                                                                                                                 
# [15] "other attached packages:"                                                                                                                                                                                                         
# [16] "[1] scales_1.0.0  tidyr_0.8.2   cowplot_0.9.3 ggplot2_3.1.0"                                                                                                                                                                      
# [17] ""                                                                                                                                                                                                                                 
# [18] "loaded via a namespace (and not attached):"                                                                                                                                                                                       
# [19] " [1] Rcpp_1.0.0       withr_2.1.2      crayon_1.3.4     dplyr_0.7.8      assertthat_0.2.0 grid_3.5.1       plyr_1.8.4       R6_2.3.0         gtable_0.2.0     magrittr_1.5     pillar_1.3.0     rlang_0.3.0.1    lazyeval_0.2.1  "
# [20] "[14] bindrcpp_0.2.2   labeling_0.3     glue_1.3.0       purrr_0.2.5      munsell_0.5.0    compiler_3.5.1   pkgconfig_2.0.2  colorspace_1.3-2 tidyselect_0.2.5 bindr_0.1.1      tibble_1.4.2    "             


#data

dat1 <- read.csv("all_methylKit_fix_gene_level_glm_tidied_with_locID_with_posinfo.csv", header = T)

## if using the expected_output file (rather than that produced by running the code) use  
## dat1 <- read.csv("./Data/expected_output/all_methylKit_fix_gene_level_glm_tidied_with_locID_with_posinfo.csv", header = T)

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




#############################################################################################
# plot along chromosome

treat_df <- as.data.frame(cbind(
c(dat2$treat5aza_dC_24h,dat2$treat5aza_dC_48h),
c(as.character(dat2$gene_name),as.character(dat2$gene_name)),
c(rep("treat5aza_dC_24h", length(dat2$treat5aza_dC_24h)), rep("treat5aza_dC_48h", length(dat2$treat5aza_dC_48h))),
c(dat2$chrID,dat2$chrID),
c(dat2$pos,dat2$pos)
))

head(treat_df )
colnames(treat_df) <- c("coeff", "gene_name", "treat", "chrID", "pos")
treat_df$coeff <-  as.numeric(as.character(treat_df$coeff))
treat_df$pos <-  as.numeric(as.character(treat_df$pos))
str(treat_df)

treat_df_chr1 <- subset(treat_df, treat_df$chrID == "1")
treat_df_chr2 <- subset(treat_df, treat_df$chrID == "2")
treat_df_chr3 <- subset(treat_df, treat_df$chrID == "3")
treat_df_chr4 <- subset(treat_df, treat_df$chrID == "4")
treat_df_chr5 <- subset(treat_df, treat_df$chrID == "5")

head(treat_df)
head(treat_df_chr1)
head(treat_df_sig_treat)

max(treat_df$coeff)
min(treat_df$coeff)

chr1_P1 <- ggplot(data=treat_df_chr1,aes(x=pos/1000000, y=coeff,colour=treat)) + 
           geom_line(aes(colour=treat, size=0.2), size=0.2) + 
           geom_hline(yintercept=0, size=0.2) + 
           ggtitle("Chromosome 1")  + xlab("Position") + ylab("log odds of methylation")  + theme(legend.position="none") + ylim(c(-24,24)) +
           scale_colour_manual(values=c("#56B4E9", "#E69F00")) +  scale_x_continuous(name="Position (Mbp)", labels = comma)

 
           
chr2_P1 <- ggplot(data=treat_df_chr2,aes(x=pos/1000000, y=coeff,colour=treat)) + 
           geom_line(aes(colour=treat), size=0.2) + 
           geom_hline(yintercept=0, size=0.2) + 
           ggtitle("Chromosome 2")  + xlab("Position") + ylab("log odds of methylation")  + theme(legend.position="none") + ylim(c(-24,24)) +
           scale_colour_manual(values=c("#56B4E9", "#E69F00"))  +  scale_x_continuous(name="Position (Mbp)", labels = comma)

            
chr3_P1 <- ggplot(data=treat_df_chr3,aes(x=pos/1000000, y=coeff,colour=treat)) + 
           geom_line(aes(colour=treat), size=0.2) + 
           geom_hline(yintercept=0, size=0.2) + 
           ggtitle("Chromosome 3")  + xlab("Position") + ylab("log odds of methylation")   + theme(legend.position="none") + ylim(c(-24,24))+
           scale_colour_manual(values=c("#56B4E9", "#E69F00"))  +  scale_x_continuous(name="Position (Mbp)", labels = comma)

           
chr4_P1 <- ggplot(data=treat_df_chr4,aes(x=pos/1000000, y=coeff,colour=treat)) + 
           geom_line(aes(colour=treat), size=0.2) + 
           geom_hline(yintercept=0, size=0.2) + 
           ggtitle("Chromosome 4")  + xlab("Position") + ylab("log odds of methylation")   + theme(legend.position="none") + ylim(c(-24,24)) +
           scale_colour_manual(values=c("#56B4E9", "#E69F00"))  +  scale_x_continuous(name="Position (Mbp)", labels = comma)


chr5_P1 <- ggplot(data=treat_df_chr5,aes(x=pos/1000000, y=coeff,colour=treat)) + 
           geom_line(aes(colour=treat), size=0.2) + 
           geom_hline(yintercept=0, size=0.2) + 
           ggtitle("Chromosome 5")  + xlab("Position") + ylab("log odds of methylation")   + theme(legend.position="none") + ylim(c(-24,24)) +
           scale_colour_manual(values=c("#56B4E9", "#E69F00"))  +  scale_x_continuous(name="Position (Mbp)", labels = comma)


pdf("pos_chr_2.pdf", width = 10, height = 15)
plot_grid(chr1_P1, chr2_P1, chr3_P1, chr4_P1, chr5_P1, ncol = 1, nrow = 5)
dev.off()
getwd() ## where has my plot gone....



########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "BSseq_analysis_of_glm_gene_level.R_sessionInfo.txt")

