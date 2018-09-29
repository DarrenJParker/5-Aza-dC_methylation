### add_NB_gene_info.py

import sys
import os
import getopt

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:t:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = None
trans_filename = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** add_NB_gene_info.py | Written by DJP, 10/08/18 in Python 3.5 in Edinburgh ****\n")
		print("Adds Refseq info to all_methylKit_fix_gene_level_glm_tidied.csv. See 1_data_prep.sh for info")
		
		print("***** USAGE *****\n")		
		print("python3 add_NB_gene_info.py -t Data/Nvit_OGSv1.2_official_id_map.txt -i all_methylKit_fix_gene_level_glm_tidied.csv")
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-t'):
		trans_filename  = arg
	else:
		print("i dont know")
		sys.exit(2)

####### get trans info into dict

trans_dict = {}

trans_file = open(trans_filename)
for line in trans_file:
	line = line.rstrip("\n").split("\t")
	gene_name = line[1]
	loc_ID = line[5].split(":")
	if len(loc_ID) == 2:
		loc_ID = loc_ID[1]
		#print(gene_name + "\t" + loc_ID)
		trans_dict[gene_name] = loc_ID


### match up (and tidy in file)

out_file = open(in_file_name.rstrip(".csv") + "_with_locID.csv" , "w")
line_N = 0
in_file = open(in_file_name)
for line in in_file:
	line = line.rstrip("\n").split(",")
	line_N = line_N + 1
	if line_N == 1:
		header = "gene_name,N_loc,Refseq_LocID"
		for i in range(3,len(line)):
			header = header + "," + line[i]
		out_file.write(header + "\n")
	else:
		gene_want = line[1]
		loc_ID = trans_dict.get(gene_want)
		if loc_ID == None:
			loc_ID = "NA"
		new_line = gene_want + "," + line[2] + "," + loc_ID
		
		for i in range(3,len(line)):
			new_line = new_line + "," + line[i]
		
		out_file.write(new_line + "\n")
		
		# print(new_line)
		# print(line)
		# 
		
		
	#print(line)
	
print("\nFinished Jonathan Harker\n\n")
