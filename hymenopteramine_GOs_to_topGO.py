## hymenopteramine_GOs_to_topGO.py

import sys
import os
import getopt


try:
	opts, args = getopt.getopt(sys.argv[1:], 'l:g:d:o:h')
																					 	
except getopt.GetoptError:
	print('ERROR getting options, please see help by specifing -h')
	sys.exit(2) ### close ofter error from try
	
	
arg_len = len(opts)

if arg_len == 0:
	print('No options provided. Please see help by specifing -h')
	sys.exit(2)


in_list_name = None
GO_term_file_name = None
DROSO_GO_term_file_name = None
ortho_file_name = None

#print (opts) ## see args
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print("\n**** hymenopteramine_GOs_to_topGO.py | Written by DJP 27/09/18 in Python 3.5 in Lausanne, Swiss ****\n")
		print("Takes datafiles from hymenopteramine and converts them for use in TopGO")
		print("\n**** USAGE ****\n")
		print("python /Users/dparker/Documents/University/St_Andrews_work/Nicki_BSseq/Nasonia_BSseq_code/hymenopteramine_GOs_to_topGO.py \ ")
		print("-g Nvit_1.2_GOs_from_hymenopteramine_270918.tsv \ ")
		print("-l ../gff/Nasonia_vitripennis.Nvit_2.1_genenames.txt \ ")
		print("-d Drosp_GOs_from_hymenopteramine_for_Nvit_1.2_orths_270918.tsv \ ")
		print("-o Gene_Orthologues_for_Nvit_1.2_from_hymenopteramine_270918.tsv\n\n")

		sys.exit(2)
		
	elif opt in ('-l'):
		in_list_name = arg
	elif opt in ('-g'):
		GO_term_file_name = arg
	elif opt in ('-d'):
		DROSO_GO_term_file_name = arg
	elif opt in ('-o'):
		ortho_file_name = arg
	else:
		print("i dont know")
		sys.exit(2)


##### get all go_term info into a dict

GO_term_file = open(GO_term_file_name)

GO_term_dict = {}

seen_gene = set()
for line in GO_term_file:
	line = line.rstrip("\n").split("\t")
	gene_ID = line[0]
	GO_term = line[1]
	
	if gene_ID not in seen_gene:
		seen_gene.add(gene_ID)
		GO_term_dict[gene_ID] = GO_term
		
	else:
		rec = GO_term_dict.get(gene_ID)
		rec = rec + ", " + GO_term
		GO_term_dict[gene_ID] = rec

# print(GO_term_dict)

##### output

in_list = open(in_list_name)
out_file = open(GO_term_file_name + "_fortopgo.txt", "w")

for line in in_list:
	curr_id = line.rstrip("\n")
	curr_GO_term_rec = GO_term_dict.get(curr_id)
	if curr_GO_term_rec == None:
		curr_GO_term_rec = ""
	
	out_file.write(curr_id + "\t" + curr_GO_term_rec + "\n")

in_list.close()
####### also do Droso?

if DROSO_GO_term_file_name != None:
	
	
	## get orth info into dict
	
	orth_info_dict = {}
	N_genes_with_mel_orths = set()

	
	seen_gene = set()
	ortho_file = open(ortho_file_name)
	for line in ortho_file:
		line = line.rstrip("\n").split("\t")
		#print(line)
		if line[3] == "D. melanogaster":

			if line[0] not in seen_gene:
				seen_gene.add(line[0])
				orth_info_dict[line[0]] = [line[2]]
			else:
				rec = orth_info_dict.get(line[0])
				rec.append(line[2])
				orth_info_dict[line[0]] = rec


	##### get all go_term info into a dict
	
	DROSO_GO_term_file = open(DROSO_GO_term_file_name)
	
	DROSO_GO_term_dict = {}
	
	seen_gene = set()
	for line in DROSO_GO_term_file :
		line = line.rstrip("\n").split("\t")
		gene_ID = line[0]
		GO_term = line[1]
		
		if gene_ID not in seen_gene:
			seen_gene.add(gene_ID)
			DROSO_GO_term_dict[gene_ID] = [GO_term]
			
		else:
			rec = DROSO_GO_term_dict.get(gene_ID)
			rec.append(GO_term)
			DROSO_GO_term_dict[gene_ID] = rec
			
	

	in_list = open(in_list_name)
	out_file_dros_inter = open(DROSO_GO_term_file_name + "_fortopgo_intersection.txt", "w")
	out_file_dros_union = open(DROSO_GO_term_file_name + "_fortopgo_union.txt", "w")	
	for line in in_list:
		curr_id = line.rstrip("\n")
		orth_list = orth_info_dict.get(curr_id)
		
		intersection_GO_list = []
		union_GO_list = set()
		
		if orth_list == None:
			intersection_GO_list = []
			union_GO_list = set()
		elif len(orth_list) == 1:
			GO_T = DROSO_GO_term_dict.get(orth_list[0])
			if GO_T == None:
				intersection_GO_list = []
				union_GO_list = set()			
			else:
				intersection_GO_list = GO_T
				union_GO_list    = set(GO_T)
		else:
			GO_term_lists = []
			for el in orth_list:
				GO_terms = DROSO_GO_term_dict.get(el)
				if GO_terms != None:
					GO_term_lists.append(GO_terms)

			if len(GO_term_lists) == 0:
				intersection_GO_list = []
				union_GO_list = set()
			elif len(GO_term_lists) >= 1:
				# print("lll")
				# print(GO_term_lists)
				GO_inter = set(GO_term_lists[0])
				for s in GO_term_lists[1:]:
					GO_inter.intersection_update(s)
				intersection_GO_list = GO_inter
				
				full_GO_list = []
				for l in GO_term_lists:
					for e in l:
						full_GO_list.append(e)
				
				union_GO_list = set(full_GO_list)	
		
		intersection_GO_list_out = ""
		for go in list(intersection_GO_list):
			intersection_GO_list_out = intersection_GO_list_out + ", " + go
		
		union_GO_list_out = ""
		for go in list(union_GO_list):
			union_GO_list_out = union_GO_list_out + ", " + go
		
		intersection_GO_list_out = intersection_GO_list_out.lstrip(", ")		
		union_GO_list_out = union_GO_list_out.lstrip(", ")	
		
		out_file_dros_inter.write(curr_id + "\t" + intersection_GO_list_out + "\n")
		out_file_dros_union.write(curr_id + "\t" + union_GO_list_out + "\n")
		
print("\n\nFinished, Joseph Banks\n\n")






