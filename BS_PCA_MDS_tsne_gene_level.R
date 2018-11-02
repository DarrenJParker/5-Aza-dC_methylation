## libs

library(ggplot2)
library(Rtsne)
library(tidyr
library(ggplot2)
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
# [16] "[1] cowplot_0.8.0 Rtsne_0.13    ggplot2_2.2.1"                                                                                                                                                  
# [17] ""                                                                                                                                                                                               
# [18] "loaded via a namespace (and not attached):"                                                                                                                                                     
# [19] " [1] colorspace_1.3-2 scales_0.5.0     compiler_3.4.1   lazyeval_0.2.1   plyr_1.8.4       gtable_0.2.0     tibble_1.3.4     Rcpp_0.12.13     grid_3.4.1       rlang_0.1.6      munsell_0.4.3   "


######################################################################
## DATA
######################################################################

input_file_name = "all_methylKit_fix_gene_level_AMFT_props.csv" #### specify name of input file here
dat1 = read.csv(input_file_name, check.names=FALSE, stringsAsFactors=FALSE) 
dat_sample_info = read.csv("Data/BSseq_sample_info.csv", check.names=FALSE, stringsAsFactors=FALSE)

### ditch NAs (removes loci with a cov of 0 in any sample)
dat2 = na.omit(dat1)

## total loci

sum(dat2$Nloci) 

## 2087098

### transpose data
dat2_t <- as.data.frame(t(dat2[3:length(colnames(dat2))]))

## PCA

PCA_res2 <- prcomp(dat2_t) 
summary(PCA_res2)
plot(PCA_res2, type = "l")

######################################################################
## PCA
######################################################################


PCA_res2 <- prcomp(dat2_t) 
summary(PCA_res2)

#### PLOT 

plot_PCA_1 <- function(PCA_output,sample_inf,tit_inf){
	PCs2 = as.data.frame(PCA_output$x)
	PCs2 = cbind(PCs2,sample_inf)
	P1 <- ggplot(PCs2, aes(x=PC1, PC2, col=treat, shape=tiss_type, label = sample_ord)) + geom_point(size = 4) + 
		theme_bw() +
		scale_colour_manual(values=c("#56B4E9", "#E69F00",  "#000000"))	+ 
		ggtitle(paste("PCA, ", tit_inf))	
	return(P1)
	
}

plot_PCA_2 <- function(PCA_output,sample_inf,tit_inf){
	PCs2 = as.data.frame(PCA_output$x)
	PCs2 = cbind(PCs2,sample_inf)
	P1 <- ggplot(PCs2, aes(x=PC1, PC2, col=treat, shape=collection, label = sample_ord)) + geom_point(size = 4) + 
		theme_bw() +
		scale_colour_manual(values=c("#56B4E9", "#E69F00",  "#000000"))	+ 
		ggtitle(paste("PCA, ", tit_inf))	
	return(P1)
	
}

P_PCA_v1   <- plot_PCA_1(PCA_res2,dat_sample_info,  paste("N Genes =", length(dat2[,1])))
P_PCA_v2   <- plot_PCA_2(PCA_res2,dat_sample_info,  paste("N Genes =", length(dat2[,1])))




######################################################################
## T-SNE
######################################################################


data_tsne = dat2_t 

set.seed(42) ### for consistency with the paper version
tsne_model_1 = Rtsne(as.matrix(data_tsne), check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)


### PLOT

plot_tsne_1 <- function(tsne_mod, sample_inf, tit_inf){
	
	## getting the two dimension matrix
	d_tsne = as.data.frame(tsne_mod$Y)

	df_tsne = as.data.frame(cbind(d_tsne, sample_inf))
	P1 <- ggplot(df_tsne, aes(x=V1, y = V2, col=treat, shape=tiss_type, label = sample_ord)) + geom_point(size = 4) + 
		theme_bw() +
		scale_colour_manual(values=c("#56B4E9", "#E69F00",  "#000000"))	+ 
		ggtitle(paste("t-SNE, ", tit_inf))	 + 
		xlab("Arbitrary units") + 
		ylab("Arbitrary units")
	
	return(P1)
}


plot_tsne_2 <- function(tsne_mod, sample_inf, tit_inf){
	
	## getting the two dimension matrix
	d_tsne = as.data.frame(tsne_mod$Y)

	df_tsne = as.data.frame(cbind(d_tsne, sample_inf))
	P1 <- ggplot(df_tsne, aes(x=V1, y = V2, col=treat, shape=collection, label = sample_ord)) + geom_point(size = 4) + 
		theme_bw() +
		scale_colour_manual(values=c("#56B4E9", "#E69F00",  "#000000"))	+ 
		ggtitle(paste("t-SNE, ", tit_inf))	 +
		xlab("Arbitrary units") + 
		ylab("Arbitrary units")	
	return(P1)
}


P_tsne_alls_v1   <- plot_tsne_1(tsne_model_1, dat_sample_info, paste("N Genes =", length(dat2[,1])))
P_tsne_alls_v2   <- plot_tsne_2(tsne_model_1, dat_sample_info, paste("N Genes =", length(dat2[,1])))


######## output

dir.create("Gene_level_plots")
setwd("Gene_level_plots")

pdf("PCA_Gene_level.pdf", width = 21, height = 10.5)
plot_grid(P_PCA_v1, P_PCA_v2, labels=c("A", "B"), ncol = 2, nrow = 1)
dev.off()
getwd() ## where has my plot gone....?

pdf("tSNE_Gene_level.pdf", width = 21, height = 10.5)
plot_grid(P_tsne_alls_v1 , P_tsne_alls_v2, labels=c("A", "B"), ncol = 2, nrow = 1)
dev.off()
getwd() ## where has my plot gone....?






########################################################################################################################################################################
####### output session info

writeLines(capture.output(sessionInfo()), "BS_PCA_MDS_tsne_gene_level_.R_sessionInfo.txt")




