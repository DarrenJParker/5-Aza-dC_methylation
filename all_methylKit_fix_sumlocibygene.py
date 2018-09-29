##### all_methylKit_fix_tidier.py

import sys
import os
import getopt
import decimal
import numpy as np

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:a:o:d:f:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = "NOTHINGSET"
annot_file_name = "NOTHINGSET"
out_prefix = "NOTHINGSET"
max_dist = 1000
cov_filt = None
#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** all_methylKit_fix_sumlocibygene.py | Written by DJP, 21/07/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Takes all_methylKit_fix.csv and annotation file (all_methylKit_fix_gene_closestgene_exon_CDS_5UTR_3UTR_annot.txt) and sums loci counts up by gene. \nAlso outputs the loci kept in a sep file.")
		
		print("\n***** USAGE *****\n")		
		print("\npython all_methylKit_fix_sumlocibygene.py -i [path to all_methylKit_fix.csv] -o [output prefix] -a [path to annotation file] [options]\n")
		print("\n***** OPTIONS *****\n")			
		print("\n-d\t[max distance] - distance (in bp) away from gene for a loci to be counted. I.e -d 100 sums counts from all loci 100 bp up- and down- stream of the gene annotation. Default: 1000")	
		print("\n-f\t[cov filter] - filter loci with <x coverage BEFORE summing up. Default: 0\n\n")	
				
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-a'):
		annot_file_name = arg
	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-d'):
		max_dist = arg
	elif opt in ('-f'):
		cov_filt = arg
	else:
		print("i dont know")
		sys.exit(2)

if cov_filt == None:
	print("\nNo cov filter set, using all SNPs\n")
else:
	cov_filt = int(cov_filt)
	print("\ncov filter set to exclude loci with <"+ str(cov_filt) + " in any sample\n")



####### get annot info + list of SNPs I want to use (this is to cut down on mem usage, as it is likely this will be a minority of SNPs in the file..)

try:
	max_dist = int(max_dist)
except:
	print("\n\nmax_dist value set to " + str(max_dist) + " which is not an integer. Please fix. Exiting\n\n\n" )
	sys.exit(2)


print("\n\nmax distance set to " + str(max_dist) + ". See help (-h) for details")

annot_file = open(annot_file_name)

kept_SNPs = 0
gene_to_SNPs_dict = {}
seen_gene = set()
seen_SNP = set()

for line in annot_file:
	line = line.rstrip("\n").split("\t")
	abs_dist = line[5]
	gene_overlap_stat = line[2]
	SNP_name = line[0]
	gene_name = line[3]
	if abs_dist != "NA":
		abs_dist = int(abs_dist)
		if abs_dist <= max_dist:
			if gene_overlap_stat == "NA": ### filter SNPs that come from overlpping gene annots
				kept_SNPs = kept_SNPs + 1
				seen_SNP.add(SNP_name)
				if gene_name not in seen_gene:
					seen_gene.add(gene_name)
					gene_to_SNPs_dict[gene_name] = set([SNP_name])
				else:
					snp_rec = gene_to_SNPs_dict.get(gene_name)
					snp_rec.add(SNP_name)
					gene_to_SNPs_dict[gene_name] = snp_rec

				
				
				
print("\nNumber of Loci kept: " + str(kept_SNPs))
print("\nNumber of genes with at least 1 selected Loci in; " + str(len(seen_gene)))


###### add SNP info into dict

in_file = open(in_file_name)
line_N = 0
sample_order = ""
SNP_dict = {}

# cov by sample

samp_cov_dict = {}
seen_sample = set()

for line in in_file:
	line_N = line_N + 1
	line = line.rstrip("\n").split(",")
	
	## get sample names in order
	if line_N == 1:
		for i in range(2, len(line)):
			sample_order = sample_order + "," + line[i]

	else:
		
		locus_name = "LN__" + line[0] + "-" + line[1]
		#print(locus_name)	
		if locus_name in seen_SNP:
			new_line = []
			for i in range(2, len(line)):
				vals = line[i].split("/")
				
				### check all count values are integer - if not throw error and exit
				
				try:
					meth_c   = int(vals[0])
					unmeth_c = int(vals[1])
				except:
					print("Not all count values are integer - exiting!")
					sys.exit(2)
				
				new_line.append(meth_c)
				new_line.append(unmeth_c)
			SNP_dict[locus_name] = new_line

in_file.close()


## filter snps with low cov out (do not add them)

SNPs_with_too_low_cov = set()

if cov_filt != None:

	for el in gene_to_SNPs_dict:
		SNP_list = gene_to_SNPs_dict.get(el)
		for s in SNP_list:
			count_list = SNP_dict.get(s)
			#print(count_list)
			
			for i in range(0,len(count_list),2):
				SNP_cov = count_list[i] + count_list[i+1]
				#print(SNP_cov)
				if SNP_cov < cov_filt:
					SNPs_with_too_low_cov.add(s)
			
print("Number of SNPs exculded as have < " + str(cov_filt) + " cov in any sample = " + str(len(SNPs_with_too_low_cov)))

###### sum SNP counts by gene

summed_counts_by_gene_dict = {}
N_snps_in_gene = {}

seen_el = set()
for el in gene_to_SNPs_dict:
	SNP_list = gene_to_SNPs_dict.get(el)
	
	N_gene_snps = 0
	for s in SNP_list:
		count_list = SNP_dict.get(s)
		# print(el)
		# print(s)
		#print(count_list)
		
		if s not in SNPs_with_too_low_cov:
			N_gene_snps = N_gene_snps + 1
			if el not in seen_el:
				seen_el.add(el)
				summed_counts_by_gene_dict[el] = count_list
			else:
				curr_snp_count = summed_counts_by_gene_dict.get(el)
				updated_snp_count = np.add(count_list, curr_snp_count).tolist()
				# print("UP")
				# print(updated_snp_count)
				summed_counts_by_gene_dict[el] = updated_snp_count
	N_snps_in_gene[el] = N_gene_snps	
	
####### output in same style as all_methylKit_fix.csv

out_file = open(out_prefix + "_gene_level.csv", "w")

header = "genename,Nloci" + sample_order

out_file.write(header + "\n")

for g in summed_counts_by_gene_dict:
	count_list = summed_counts_by_gene_dict.get(g)
	curr_nth = 0
	new_line = ""
	SNP_N = N_snps_in_gene.get(g)

	for c in count_list:
		curr_nth = curr_nth + 1
		if curr_nth % 2 == 0:
			new_line = new_line + "/" + str(c)
		else:
			new_line = new_line + "," + str(c)
	#print(g + "," + str(SNP_N) + new_line)
	out_file.write(g + "," + str(SNP_N) + new_line + "\n")


######## output SNPs in the summed gene file


out_file_SNPs_in_genes = open(out_prefix + "_Loci_in_genes.csv", "w")
SNPs_in_genes_kept = 0
in_file = open(in_file_name)
line_N = 0
for line in in_file:
	line_N = line_N + 1
	line1 = line.rstrip("\n")
	line = line.rstrip("\n").split(",")
	
	## get sample names in order
	if line_N == 1:
		out_file_SNPs_in_genes.write(line1 + "\n")

	else:
		locus_name = "LN__" + line[0] + "-" + line[1]
		
		
		if locus_name in seen_SNP:
			if locus_name not in SNPs_with_too_low_cov:
				out_file_SNPs_in_genes.write(line1 + "\n")
				SNPs_in_genes_kept = SNPs_in_genes_kept + 1
in_file.close()


print("Number of Loci in genes kept = " + str(SNPs_in_genes_kept) + "\n")



print("Finished, Jimmy Parker\n\n\n")


