#####  BSseq_binomial_glm_tidier.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'd:o:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_dir_name = "NOTHINGSET"
out_prefix = "NOTHINGSET"


#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n****  BSseq_binomial_glm_tidier.py | Written by DJP, 17/06/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Running BSseq_binomial_glm.R gives two files per locus: a _LLChisq.csv and _coeff.csv file")
		
		sys.exit(2)
		
	elif opt in ('-d'):
		in_dir_name = arg
	elif opt in ('-o'):
		out_prefix = arg	
	else:
		print("i dont know")
		sys.exit(2)



###### read all files from input dir join together in a useful format

## annoyingly when the coeffs are exported - not all are as if there were not enough for a partic group it is not outputted ...
## Going to add all info I do have to sep dicts - then pull this putting NA when == None
## This also means I have to specify the coeffs beforehand (I could run through all files - but there are > 4,000,000...)

# model coeffs I want:
# 	
# (Intercept)
# dat1$treat5aza_dC_24h
# dat1$treat5aza_dC_48h
# dat1$collection24h
# dat1$collection48h
# dat1$tiss_typeHeads
# dat1$treat5aza_dC_24h:dat1$collection24h
# dat1$treat5aza_dC_48h:dat1$collection24h
# dat1$treat5aza_dC_24h:dat1$collection48h
# dat1$treat5aza_dC_48h:dat1$collection48h
# dat1$treat5aza_dC_24h:dat1$tiss_typeHeads
# dat1$treat5aza_dC_48h:dat1$tiss_typeHeads
# dat1$collection24h:dat1$tiss_typeHeads
# dat1$collection48h:dat1$tiss_typeHeads
# dat1$treat5aza_dC_24h:dat1$collection24h:dat1$tiss_typeHeads
# dat1$treat5aza_dC_48h:dat1$collection24h:dat1$tiss_typeHeads
# dat1$treat5aza_dC_24h:dat1$collection48h:dat1$tiss_typeHeads
# dat1$treat5aza_dC_48h:dat1$collection48h:dat1$tiss_typeHeads

## stats i want
# dat1$treat
# dat1$collection
# dat1$tiss_type
# dat1$treat:dat1$collection
# dat1$treat:dat1$tiss_type
# dat1$collection:dat1$tiss_type
# dat1$treat:dat1$collection:dat1$tiss_type
# 

### coeff dicts (log-odd ratios)

intercept_dict = {}

treat5aza_dC_24h_dict = {}
treat5aza_dC_48h_dict = {}

collection24h_dict = {}
collection48h_dict = {}

tiss_typeHeads_dict = {}

treat5aza_dC_24h_collection24h_dict = {}
treat5aza_dC_48h_collection24h_dict = {}
treat5aza_dC_24h_collection48h_dict = {}
treat5aza_dC_48h_collection48h_dict = {}
									
treat5aza_dC_24h_tiss_typeHeads_dict = {}
treat5aza_dC_48h_tiss_typeHeads_dict = {}
									
collection24h_tiss_typeHeads_dict = {}
collection48h_tiss_typeHeads_dict = {}
									
treat5aza_dC_24h_collection24h_tiss_typeHeads_dict = {}
treat5aza_dC_48h_collection24h_tiss_typeHeads_dict = {}
treat5aza_dC_24h_collection48h_tiss_typeHeads_dict = {}
treat5aza_dC_48h_collection48h_tiss_typeHeads_dict = {}

## stat dicts

LL_treat_dict = {}
LL_collection_dict = {}
LL_tiss_type_dict = {}
LL_treat_collection_dict = {}
LL_treat_tiss_type_dict = {}
LL_collection_tiss_type_dict = {}
LL_treat_collection_tiss_type_dict = {}


locus_l = []

path = in_dir_name
for path, subdirs, files in os.walk(path):
	for name in files:
		if name.endswith("_coeff.csv"):
			#print (os.path.join(path, name))
			curr_coeff_f = open(os.path.join(path, name))
			
			locus_name = name.split(".csv")[0]
			#print(locus_name )
			locus_l.append(locus_name )
			for line in curr_coeff_f:
				line = line.rstrip("\n").split(",")
				coeff_name = line[0].strip('"')
				coeff_est  = line[1].strip('"')
				
				if coeff_name == "(Intercept)":
					intercept_dict[locus_name] = coeff_est
				
				if coeff_name == "dat1$treat5aza_dC_24h":
					treat5aza_dC_24h_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_48h":
					treat5aza_dC_48h_dict[locus_name] = coeff_est

				if coeff_name == "dat1$collection24h":
					collection24h_dict[locus_name] = coeff_est
				if coeff_name == "dat1$collection48h":
					collection48h_dict[locus_name] = coeff_est
					
				if coeff_name == "dat1$tiss_typeHeads":
					tiss_typeHeads_dict[locus_name] = coeff_est
				
				if coeff_name == "dat1$treat5aza_dC_24h:dat1$collection24h":
					treat5aza_dC_24h_collection24h_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_48h:dat1$collection24h":
					treat5aza_dC_48h_collection24h_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_24h:dat1$collection48h":
					treat5aza_dC_24h_collection48h_dict[locus_name] = coeff_est				
				if coeff_name == "dat1$treat5aza_dC_48h:dat1$collection48h":
					treat5aza_dC_48h_collection48h_dict[locus_name] = coeff_est

				if coeff_name == "dat1$treat5aza_dC_24h:dat1$tiss_typeHeads":
					treat5aza_dC_24h_tiss_typeHeads_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_48h:dat1$tiss_typeHeads":
					treat5aza_dC_48h_tiss_typeHeads_dict[locus_name] = coeff_est				
					
				if coeff_name == "dat1$collection24h:dat1$tiss_typeHeads":
					collection24h_tiss_typeHeads_dict[locus_name] = coeff_est				
				if coeff_name == "dat1$collection48h:dat1$tiss_typeHeads":
					collection48h_tiss_typeHeads_dict[locus_name] = coeff_est					
					
				if coeff_name == "dat1$treat5aza_dC_24h:dat1$collection24h:dat1$tiss_typeHeads":
					treat5aza_dC_24h_collection24h_tiss_typeHeads_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_48h:dat1$collection24h:dat1$tiss_typeHeads":
					treat5aza_dC_48h_collection24h_tiss_typeHeads_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_24h:dat1$collection48h:dat1$tiss_typeHeads":
					treat5aza_dC_24h_collection48h_tiss_typeHeads_dict[locus_name] = coeff_est
				if coeff_name == "dat1$treat5aza_dC_48h:dat1$collection48h:dat1$tiss_typeHeads":
					treat5aza_dC_48h_collection48h_tiss_typeHeads_dict[locus_name] = coeff_est
			
			curr_coeff_f.close()
		
			curr_LL_f = open(os.path.join(path, name.replace("_coeff.csv", "_LLChisq.csv")))
			#print(os.path.join(path, name.replace("_coeff.csv", "_LLChisq.csv")))
			for line in curr_LL_f:
				line = line.rstrip("\n").split(",")
				term_name = line[0].strip('"')
				chi_val   = line[4].strip('"')
				pval      = line[5].strip('"')
				
				if term_name == "dat1$treat":
					LL_treat_dict[locus_name] = chi_val + "," + pval
				if term_name == "dat1$collection":
					LL_collection_dict[locus_name] = chi_val + "," + pval
				if term_name == "dat1$tiss_type":
					LL_tiss_type_dict[locus_name] = chi_val + "," + pval
				if term_name == "dat1$treat:dat1$collection":
					LL_treat_collection_dict[locus_name] = chi_val + "," + pval
				if term_name == "dat1$treat:dat1$tiss_type":
					LL_treat_tiss_type_dict[locus_name] = chi_val + "," + pval
				if term_name == "dat1$collection:dat1$tiss_type":
					LL_collection_tiss_type_dict[locus_name] = chi_val + "," + pval
				if term_name == "dat1$treat:dat1$collection:dat1$tiss_type":
					LL_treat_collection_tiss_type_dict[locus_name] = chi_val + "," + pval

			curr_LL_f.close()
				

#### sort the loci (just for tidiness...) --> then export

output_file_name = out_prefix + "_glm_tidied.csv" 
output_file = open(output_file_name, "w")
output_file.write(
		"Locus_name" + "," + "chrom" + "," + "pos" + "," +
		"intercept" +  "," +
		"treat5aza_dC_24h" + "," + "treat5aza_dC_48h" + "," + "LRT_treat,p_treat" + "," +	
		"collection24h" + "," + "collection48h" + "," + "LRT_collection,p_collection" + ","  +
		"tiss_typeHeads" + "," + "LRT_tiss_type,p_tiss_type" + "," +
		"treat5aza_dC_24h_collection24h" + "," + "treat5aza_dC_48h_collection24h" + "," + "treat5aza_dC_24h_collection48h" + "," + "treat5aza_dC_48h_collection48h" + "," + "LRT_treat_by_collection,p_treat_by_collection" + "," +
		"treat5aza_dC_24h_tiss_typeHeads" + "," + "treat5aza_dC_48h_tiss_typeHeads"  + "," +  "LRT_treat_by_tiss_type,p_treat_by_tiss_type" + "," +
		"collection24h_tiss_typeHeads" + "," + "collection48h_tiss_typeHeads" + "," + "LRT_collection_by_tiss_type,p_collection_by_tiss_type" + "," +
		"treat5aza_dC_24h_collection24h_tiss_typeHeads" + "," + "treat5aza_dC_48h_collection24h_tiss_typeHeads" + "," +
		"treat5aza_dC_24h_collection48h_tiss_typeHeads" + "," + "treat5aza_dC_48h_collection48h_tiss_typeHeads" + "," + "LRT_treat_by_collection_by_tiss_type,p_treat_by_collection_by_tiss_type" + "\n")



print(len(locus_l))

locus_l_sorted = sorted(locus_l)

for loc in locus_l_sorted:
	intercept = intercept_dict.get(loc)
	treat5aza_dC_24h = treat5aza_dC_24h_dict.get(loc)
	treat5aza_dC_48h = treat5aza_dC_48h_dict.get(loc)
	collection24h = collection24h_dict.get(loc)
	collection48h = collection48h_dict.get(loc)
	tiss_typeHeads = tiss_typeHeads_dict.get(loc)
	treat5aza_dC_24h_collection24h = treat5aza_dC_24h_collection24h_dict.get(loc)
	treat5aza_dC_48h_collection24h = treat5aza_dC_48h_collection24h_dict.get(loc)
	treat5aza_dC_24h_collection48h = treat5aza_dC_24h_collection48h_dict.get(loc)
	treat5aza_dC_48h_collection48h = treat5aza_dC_48h_collection48h_dict.get(loc)
	treat5aza_dC_24h_tiss_typeHeads = treat5aza_dC_24h_tiss_typeHeads_dict.get(loc)
	treat5aza_dC_48h_tiss_typeHeads = treat5aza_dC_48h_tiss_typeHeads_dict.get(loc)
	collection24h_tiss_typeHeads = collection24h_tiss_typeHeads_dict.get(loc)
	collection48h_tiss_typeHeads = collection48h_tiss_typeHeads_dict.get(loc)
	treat5aza_dC_24h_collection24h_tiss_typeHeads = treat5aza_dC_24h_collection24h_tiss_typeHeads_dict.get(loc)
	treat5aza_dC_48h_collection24h_tiss_typeHeads = treat5aza_dC_48h_collection24h_tiss_typeHeads_dict.get(loc)
	treat5aza_dC_24h_collection48h_tiss_typeHeads = treat5aza_dC_24h_collection48h_tiss_typeHeads_dict.get(loc)
	treat5aza_dC_48h_collection48h_tiss_typeHeads = treat5aza_dC_48h_collection48h_tiss_typeHeads_dict.get(loc)
	LL_treat = LL_treat_dict.get(loc)
	LL_collection = LL_collection_dict.get(loc)
	LL_tiss_type = LL_tiss_type_dict.get(loc)
	LL_treat_collection = LL_treat_collection_dict.get(loc)
	LL_treat_tiss_type = LL_treat_tiss_type_dict.get(loc)
	LL_collection_tiss_type = LL_collection_tiss_type_dict.get(loc)
	LL_treat_collection_tiss_type = LL_treat_collection_tiss_type_dict.get(loc)
	
	if intercept == None:
		intercept = 'NA'
	if treat5aza_dC_24h == None:
		treat5aza_dC_24h = 'NA'
	if treat5aza_dC_48h == None:
		treat5aza_dC_48h = 'NA'
	if collection24h == None:
		collection24h = 'NA'
	if collection48h == None:
		collection48h = 'NA'
	if tiss_typeHeads == None:
		tiss_typeHeads = 'NA'
	if treat5aza_dC_24h_collection24h == None:
		treat5aza_dC_24h_collection24h = 'NA'
	if treat5aza_dC_48h_collection24h == None:
		treat5aza_dC_48h_collection24h = 'NA'
	if treat5aza_dC_24h_collection48h == None:
		treat5aza_dC_24h_collection48h = 'NA'
	if treat5aza_dC_48h_collection48h == None:
		treat5aza_dC_48h_collection48h = 'NA'
	if treat5aza_dC_24h_tiss_typeHeads == None:
		treat5aza_dC_24h_tiss_typeHeads = 'NA'
	if treat5aza_dC_48h_tiss_typeHeads == None:
		treat5aza_dC_48h_tiss_typeHeads = 'NA'
	if collection24h_tiss_typeHeads == None:
		collection24h_tiss_typeHeads = 'NA'
	if collection48h_tiss_typeHeads == None:
		collection48h_tiss_typeHeads = 'NA'
	if treat5aza_dC_24h_collection24h_tiss_typeHeads == None:
		treat5aza_dC_24h_collection24h_tiss_typeHeads = 'NA'
	if treat5aza_dC_48h_collection24h_tiss_typeHeads == None:
		treat5aza_dC_48h_collection24h_tiss_typeHeads = 'NA'
	if treat5aza_dC_24h_collection48h_tiss_typeHeads == None:
		treat5aza_dC_24h_collection48h_tiss_typeHeads = 'NA'
	if treat5aza_dC_48h_collection48h_tiss_typeHeads == None:
		treat5aza_dC_48h_collection48h_tiss_typeHeads = 'NA'
	if LL_treat == None:
		LL_treat = 'NA'
	if LL_collection == None:
		LL_collection = 'NA'
	if LL_tiss_type == None:
		LL_tiss_type = 'NA'
	if LL_treat_collection == None:
		LL_treat_collection = 'NA'
	if LL_treat_tiss_type == None:
		LL_treat_tiss_type = 'NA'
	if LL_collection_tiss_type == None:
		LL_collection_tiss_type = 'NA'
	if LL_treat_collection_tiss_type == None:
		LL_treat_collection_tiss_type = 'NA'	
	
	real_loc = "LN_" + loc.split("LN")[1] 
	chrom =  loc.split("LN_")[1].split("-")[0]
	pos   =  loc.split("LN_")[1].split("-")[1]
	
	output_file.write(
		real_loc + "," + chrom + "," + pos + "," +
		intercept +  "," +
		treat5aza_dC_24h + "," + treat5aza_dC_48h + "," + LL_treat + "," +	
		collection24h + "," + collection48h + "," + LL_collection + ","  +
		tiss_typeHeads + "," + LL_tiss_type + "," +
		treat5aza_dC_24h_collection24h + "," + treat5aza_dC_48h_collection24h + "," + treat5aza_dC_24h_collection48h + "," + treat5aza_dC_48h_collection48h + "," + LL_treat_collection + "," +
		treat5aza_dC_24h_tiss_typeHeads + "," + treat5aza_dC_48h_tiss_typeHeads  + "," +  LL_treat_tiss_type + "," +
		collection24h_tiss_typeHeads + "," + collection48h_tiss_typeHeads + "," + LL_collection_tiss_type + "," +
		treat5aza_dC_24h_collection24h_tiss_typeHeads + "," + treat5aza_dC_48h_collection24h_tiss_typeHeads + "," +
		treat5aza_dC_24h_collection48h_tiss_typeHeads + "," + treat5aza_dC_48h_collection48h_tiss_typeHeads + "," + LL_treat_collection_tiss_type + "\n")


print("\n\nFinished, Jules\n\n")



