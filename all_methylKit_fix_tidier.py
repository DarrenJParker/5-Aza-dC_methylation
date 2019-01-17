##### all_methylKit_fix_tidier.py

import sys
import os
import getopt
import decimal


try:
	opts, args = getopt.getopt(sys.argv[1:], 'i:o:m:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_file_name = "NOTHINGSET"
out_prefix = "NOTHINGSET"
min_coverage = 0

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** all_methylKit_fix_tidier.py | Written by DJP, 05/06/18 in Python 3.5 in Lausanne, Swiss ****\n")
		
		print("Takes all_methylKit_fix.csv and outputs prop of methylated sites")
		print("File: all_methylKit_fix.csv contains methylation data for all samples in one file.")
		print("Each row corresponds to a CpG in the analysis")
		print("There are 53 columns the first two are:")
		print("\t1: CHR - name of chromosome CpG is on")
		print("\t2: POS - Chromosome co-ordinate of CpG")
		print("Each subsequent column corresponds to a sample and shows the number of reads in which this CpG appears to be methylated and the total no. of reads covering this CpG (in the form no. reads methylated/total reads)\n\n") 

		print("\n\n***** USAGE *****\n")		
		print("\n\npython all_methylKit_fix_tidier.py -i [path to all_methylKit_fix.csv] -o [output prefix]\n")
		
		
		sys.exit(2)
		
	elif opt in ('-i'):
		in_file_name = arg
	elif opt in ('-o'):
		out_prefix = arg
	elif opt in ('-m'):
		min_coverage = arg
	else:
		print("i dont know")
		sys.exit(2)

try:
	int(min_coverage)
except:
	print("\nmin_coverage was set to: " + str(min_coverage) + " which is not an integer. Exiting!\n")
	sys.exit(2)
print("\nMin coverage for any site is >= to " + str(min_coverage) + ". Sites with coverage less than this will be set to NA")



###### read file and calc proportions

in_file = open(in_file_name)
out_prop_file_name = out_prefix + "_AMFT_props.csv"
out_prop_file = open(out_prop_file_name, "w")

sites_with_0_cov = 0
sites_below_cov_thresh = 0

line_N = 0

for line in in_file:
	line_N = line_N + 1
	line = line.rstrip("\n")
	
	## output header as is
	if line_N == 1:
		out_prop_file.write(line + "\n")
	
	## split line -> take all cols after the first two and calc a prop (as a decimal NOT a float) -> gather back up -> output to file
	else:
		line = line.split(",")
		#print(line)
		
		meth_props_all = ""
		
		for i in range(2, len(line)):
			vals = line[i].split("/")
			
			### check all count values are integer - if not throw error and exit

			meth_cd = decimal.Decimal(-10000.777)
			unmeth_cd = decimal.Decimal(-10000.777)
			
			try:
				meth_cTEST = int(vals[0])
				total_cTEST = int(vals[1])
				meth_cd = decimal.Decimal(vals[0])
				total_cd = decimal.Decimal(vals[1])
			except:
				print("Not all count values are integer - exiting!")
				sys.exit(2)
			#print(vals)
			
			#### when there is 0 meth and unmeth - set meth_prop to NA
			
			meth_prop = decimal.Decimal(-50000.777)
			
			if total_cd == 0:
				meth_prop = "NA"
				sites_with_0_cov = sites_with_0_cov + 1
			elif total_cd < int(min_coverage):
				meth_prop = "NA"
				sites_below_cov_thresh = sites_below_cov_thresh + 1				
			else:
				meth_prop = decimal.Decimal(meth_cd/(total_cd))
			
			meth_props_all = meth_props_all + "," + str(meth_prop)
		
		### write props to file
		
		out_prop_file.write(line[0] + "," + line[1] + meth_props_all + "\n")

print("\nNumber of lines in " + in_file_name + " (incl. header): " + str(line_N) + "\n\n")

print("Prop file outputted to: " + out_prop_file_name + "\n\n")

print("Sites with 0 coverage: " + str(sites_with_0_cov))
print("Sites below coverage threshold (incl 0s): " + str(sites_with_0_cov + sites_below_cov_thresh))

