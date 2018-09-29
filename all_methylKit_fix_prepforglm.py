##### all_methylKit_fix_tidier.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:s:o:m:n:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = "NOTHINGSET"
sample_file_name = "NOTHINGSET" 
out_prefix = "NOTHINGSET"
min_coverage = 0
N_samples = 0

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** all_methylKit_fix_prepforglm.py | Written by DJP, 09/06/18 in Python 3.5 in Lausanne, Swiss ****\n")
		
		print("Takes all_methylKit_fix.csv or all_methylKit_fix_gene_level.csv and sep files for each locus or gene for use in R (binomial glm)")
		print("Also produces values of per locus/gene and per sample coverage\n")
		
		print("File: all_methylKit_fix.csv/all_methylKit_fix_gene_level.csv contains methylation data for all samples in one file.")
		print("Each row corresponds to a CpG in the analysis")
		print("There are 53 columns the first two are:")
		print("\t1: CHR - name of chromosome CpG is on")
		print("\t2: POS - Chromosome co-ordinate of CpG")
		print("Each subsequent column corresponds to a sample and shows the number of reads in which this CpG appears to be methylated and the total no. of reads covering this CpG (in the form no. reads methylated/total reads)\n\n") 

		print("***** USAGE *****\n")		
		print("\npython3 all_methylKit_fix_tidier.py -i [path to all_methylKit_fix.csv/all_methylKit_fix_gene_level.csv] -o [output prefix] -s [path to sample info file]\n")
		print("\n***** OPTIONS *****\n")			
		print("\n-i\t[path to all_methylKit_fix.csv] \n-o\t[output prefix]")		
		print("-s\t[path to sample info file] - should be Data/BSseq_sample_info.csv")
		print("\n\n***** Optional filtering options *****\n")
		
		print("-m\t[min coverage] - minimum coverage threshold for a sample in a locus")
		print("-mn\t[N_samples]   - number of samples -1 that can be below the min cov threshold and not be filtered out")

		print("\nfor example specifying -n 4 and -m 1 means any loci with 4 or more samples with 0 coverage will be filtered out\n\n\n")
		
		
		
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-s'):
		sample_file_name = arg
	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-m'):
		min_coverage = arg
	elif opt in ('-n'):
		N_samples = arg		
	else:
		print("i dont know")
		sys.exit(2)

filt_stat = ""

if min_coverage == 0 and N_samples == 0:
	filt_stat = "NO"
	print("\n\nNo filter set requested\n\n")
else:
	filt_stat = "YES"
	try:
		int(min_coverage)
	except:
		print("\nmin_coverage was set to: " + str(min_coverage) + " which is not an integer. Exiting!\n")
		sys.exit(2)
	try:
		int(N_samples)
	except:
		print("\nN_samples was set to: " + str(N_samples) + " which is not an integer. Exiting!\n")
		sys.exit(2)		
	
	if N_samples == 0:
		print("\n\nPlease specify how many samples the filter should be applyed to. Exiting!\n")
		sys.exit(2)		
	if min_coverage == 0:
		print("\n\nPlease specify a minimum coverage for filtering. Exiting!\n")
		sys.exit(2)				
		
	print("\n\nProducing a filtered set: Filtering loci with " + str(N_samples) + " or more samples with coverage <= to " + str(min_coverage) + " reads.\n") 


### read in sample info into a dict

line_N = 0
sample_info_dict = {}
sample_file = open(sample_file_name)
samp_header = ""

for line in sample_file:
	line_N = line_N + 1
	line = line.rstrip("\n").split(",")
	
	if line_N == 1:
		for i in range(1, len(line)):
			samp_header = samp_header + "," + line[i]
					
	else:
		vals = ""
		for i in range(1, len(line)):
			vals = vals + "," + line[i]
		
		sample_info_dict[line[0]] = vals
		# print(line)
		# print(vals)


###### read file and output

### mk output dir
output_dir_name = "single_files_for_Rglm"
try:
	os.mkdir(output_dir_name)
except OSError as exc:
	print("\n********* WARNING ************\n" + output_dir_name + " already exists. files may be overwritten.")

in_file = open(in_file_name)
line_N = 0
sample_order = []
cov_by_locus_file = open(out_prefix + "_cov_by_loc.csv", "w")
cov_by_locus_file.write("locus_name,locus_cov\n")	

# cov by sample

samp_cov_dict = {}
seen_sample = set()

for line in in_file:
	line_N = line_N + 1
	line = line.rstrip("\n").split(",")
	
	## get sample names in order
	if line_N == 1:
		for i in range(2, len(line)):
			sample_order.append(line[i])

	else:
		
		locus_name = "LN_" + line[0] + "-" + line[1]
		#print(locus_name)	
		
		out_file_name =  out_prefix + locus_name + ".csv"
		out_file_name_wp = os.path.join(output_dir_name, out_file_name)	
		out_file_wp = open(out_file_name_wp, "w")
		out_file_wp.write("Locus,meth_c,unmeth_c" + samp_header + "\n")
		
		locus_cov = 0
		
		for i in range(2, len(line)):
			vals = line[i].split("/")
			
			### check all count values are integer - if not throw error and exit
			
			try:
				meth_c   = int(vals[0])
				unmeth_c = int(vals[1])
			except:
				print("Not all count values are integer - exiting!")
				sys.exit(2)
			
			locus_cov = locus_cov + meth_c + unmeth_c
			
			sample_name_c = sample_order[i -2]
			sample_info_c = sample_info_dict.get(sample_name_c)
			
			if sample_name_c not in seen_sample:
				samp_cov = meth_c + unmeth_c
				samp_cov_dict[sample_name_c] = int(samp_cov)
				seen_sample.add(sample_name_c)
			else:
				samp_cov = samp_cov_dict.get(sample_name_c) + meth_c + unmeth_c
				samp_cov_dict[sample_name_c] = samp_cov
				
			out_file_wp.write(locus_name + "," + str(meth_c) + "," + str(unmeth_c) + sample_info_c + "\n")
		
		out_file_wp.close()	
		cov_by_locus_file.write(locus_name + "," + str(locus_cov) + "\n")	
cov_by_locus_file.close()
in_file.close()
			
# print(sample_order)
# print(samp_cov_dict)

#### output cov per samp file

cov_by_sample_file = open(out_prefix + "_cov_by_sample.csv", "w")

cov_by_sample_file.write("sample_name,total_read_counts,N_loci,coverage\n")
for el in sample_order:
	cov_sample = samp_cov_dict.get(el)
	S_cov = cov_sample / (line_N -1)
	cov_by_sample_file.write(el + "," + str(cov_sample) + "," + str(line_N -1) + "," + str(S_cov)  +  "\n")
	



########################################################################################################################################
###### Filter set 

### produce a filtered set
if filt_stat == "YES":
	
	### ID genes to be filtered out
	
	to_filter = set()
	
	### mk output dir
	output_dir_name = "single_files_for_Rglm_filtered_Nsamp_" + str(N_samples) + "Mcov_" + str(min_coverage)
	try:
		os.mkdir(output_dir_name)
	except OSError as exc:
		print("\n********* WARNING ************\n" + output_dir_name + " already exists. files may be overwritten.")
	
		
	in_file = open(in_file_name)
	line_N = 0
	
	for line in in_file:
		line_N = line_N + 1
		line = line.rstrip("\n").split(",")
		
		## get sample names in order
		if line_N > 1:
			
			locus_name = "LN_" + line[0] + "-" + line[1]
			#print(locus_name)	
			
			N_loci_below_cov_thresh = 0
			
			for i in range(2, len(line)):
				vals = line[i].split("/")
				try:
					meth_c   = int(vals[0])
					unmeth_c = int(vals[1])
				except:
					print("Not all count values are integer - exiting!")
					sys.exit(2)
				samp_cov = meth_c + unmeth_c
				
				#print(samp_cov)
				
				if samp_cov <= int(min_coverage):
					#print(samp_cov)
					N_loci_below_cov_thresh = N_loci_below_cov_thresh + 1
			
			if N_loci_below_cov_thresh >= int(N_samples):
				to_filter.add(locus_name)
	
	in_file.close()
	#print(to_filter)

	total_loci = line_N -1
	N_loci_filtered = len(to_filter)

	if total_loci == N_loci_filtered:
		print("Current filtering parameters mean all loci will be filtered! Exiting!\n\n")
		sys.exit(2)

	print("\nTotal loci in " + in_file_name + " = " + str(total_loci))
	print("Number of loci filtered = " + str(N_loci_filtered))
	filt_stats_file = open(out_prefix + "_filter_stats_" + str(N_samples) + "Mcov_" + str(min_coverage) + ".txt", "w")
	filt_stats_file.write("Total loci in " + in_file_name + " = " + str(total_loci) + "\nNumber of loci filtered = " + str(N_loci_filtered))
	filt_stats_file.close()
	
	## output filtered set

	
	in_file = open(in_file_name)
	line_N = 0
	sample_order = []
	cov_by_locus_file = open(out_prefix + "_cov_by_loc_filt" + str(N_samples) + "Mcov_" + str(min_coverage) + ".csv", "w")
	cov_by_locus_file.write("locus_name,locus_cov\n")	
	N_filt = 0
	
	# cov by sample
	
	samp_cov_dict = {}
	seen_sample = set()
	
	for line in in_file:
		line_N = line_N + 1
		line = line.rstrip("\n").split(",")
		
		## get sample names in order
		if line_N == 1:
			for i in range(2, len(line)):
				sample_order.append(line[i])
	
		else:
			
			locus_name = "LN_" + line[0] + "-" + line[1]

			#print(locus_name)	
			if locus_name not in to_filter:
				#print(locus_name)	
				N_filt = N_filt + 1
				out_file_name =  out_prefix + locus_name + ".csv"
				out_file_name_wp = os.path.join(output_dir_name, out_file_name)	
				out_file_wp = open(out_file_name_wp, "w")
				out_file_wp.write("Locus,meth_c,unmeth_c" + samp_header + "\n")
				
				locus_cov = 0
				
				for i in range(2, len(line)):
					vals = line[i].split("/")
					
					### check all count values are integer - if not throw error and exit
					
					try:
						meth_c   = int(vals[0])
						unmeth_c = int(vals[1])
					except:
						print("Not all count values are integer - exiting!")
						sys.exit(2)
					
					locus_cov = locus_cov + meth_c + unmeth_c
					
					sample_name_c = sample_order[i -2]
					sample_info_c = sample_info_dict.get(sample_name_c)
					
					if sample_name_c not in seen_sample:
						samp_cov = meth_c + unmeth_c
						samp_cov_dict[sample_name_c] = int(samp_cov)
						seen_sample.add(sample_name_c)
					else:
						samp_cov = samp_cov_dict.get(sample_name_c) + meth_c + unmeth_c
						samp_cov_dict[sample_name_c] = samp_cov
						
					out_file_wp.write(locus_name + "," + str(meth_c) + "," + str(unmeth_c) + sample_info_c + "\n")
					
				out_file_wp.close()	
				cov_by_locus_file.write(locus_name + "," + str(locus_cov) + "\n")	
	cov_by_locus_file.close()
	in_file.close()
				
	# print(sample_order)
	# print(samp_cov_dict)
	
	#### output cov per samp file
	
	cov_by_sample_file = open(out_prefix + "_cov_by_sample_filt" + str(N_samples) + "Mcov_" + str(min_coverage) + ".csv", "w")
	
	cov_by_sample_file.write("sample_name,total_read_counts,N_loci,coverage\n")
	for el in sample_order:
		cov_sample = samp_cov_dict.get(el)
		S_cov = cov_sample / N_filt
		cov_by_sample_file.write(el + "," + str(cov_sample) + "," + str(N_filt) + "," + str(S_cov)  +  "\n")


print("\n\nFinished, Dr Cook\n\n\n")



