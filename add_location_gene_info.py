### add_location_gene_info.py

import sys
import os
import getopt

try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:g:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = None
gff_filename = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** add_location_gene_info.py | Written by DJP, 09/01/19 in Python 3.5 in Lausanne ****\n")
		print("Adds chromosome and position (mid-gene) info to  all_methylKit_fix_gene_level_glm_tidied_with_locID.csv. See 1_data_prep.sh for info")
		
		print("***** USAGE *****\n")		
		print("python3 add_location_gene_info.py -g ./Data/Nasonia_vitripennis.Nvit_2.1.40_gene.gff3 -i all_methylKit_fix_gene_level_glm_tidied_with_locID.csv")
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-g'):
		gff_filename  = arg
	else:
		print("i dont know")
		sys.exit(2)

####### get gff info into dict

gff_dict_pos = {}
gff_dict_chr = {}
gff_file = open(gff_filename)
for line in gff_file:
	line = line.rstrip("\n").split("\t")
	gene_name = line[8].split(";")[0].split(":")[1]
	chr_name = line[0]
	start_pos = int(line[3])
	end_pos = int(line[4])
	mid_pos = int((end_pos + start_pos) / 2)
	gff_dict_pos[gene_name] = mid_pos
	gff_dict_chr[gene_name] = chr_name

### match up (and tidy in file)

out_file = open(in_file_name.replace(".csv", "") + "_with_posinfo.csv" , "w")
line_N = 0
in_file = open(in_file_name)
for line in in_file:
	line = line.rstrip("\n")
	line_N = line_N + 1
	if line_N == 1:
		out_file.write(line + ",chrID,pos\n")
	else:
		gene_want = line.split(",")[0]
		chr_ID = gff_dict_chr.get(gene_want)
		pos_ID = gff_dict_pos.get(gene_want)
		
		out_file.write(line + "," + chr_ID + "," + str(pos_ID) + "\n")
		
		# print(new_line)
		# print(line)
		# 
		
		
	#print(line)
	
print("\nFinished Jonathan Harker\n\n")
