import pandas as pd
import numpy as np
import glob
import skbio 
from skbio import TreeNode

fa_dir='genes_out_AZ51_filt/'
tree_files=glob.glob(fa_dir + '*tree/tree.nwk')

for tree_path in tree_files:
	tree = skbio.io.read(tree_path, format="newick", into=TreeNode)
	# get gene name 
	gene_with_tree = tree_path.split('/')[1]
	tree_index = gene_with_tree.index("_tree")
	gene_name = gene_with_tree[:tree_index]
	print(gene_name)
	# empty lists 
	node_names = []
	sample_names = []
	# make list of node names 
	for node in tree.tips():
		node_names.append(node.name)
	# make list of the sample names 
	for name in node_names:
		name = name[len(gene_name)+1:]
		expanded_name = name.split("_")
		final_index = expanded_name.index("mut")
		sample_name = "_".join(expanded_name[: final_index + 1])
		if sample_name not in sample_names:
			sample_names.append(sample_name)
	# put array together 
	empty_array = np.zeros((len(sample_names), len(node_names)))
	for i in range(len(sample_names)):
		for j in range(len(node_names)):
			if sample_names[i] in node_names[j]:
				empty_array[i][j] = 1
	biom_table = pd.DataFrame(data=empty_array, columns=node_names,
				  index=sample_names)
	biom_table = biom_table.T	
	# figure out the name 
	pf_index = tree_path.find("AZ51_filt/") + 10
	tree_index = tree_path.find("_tree")
	pf_name = tree_path[pf_index:tree_index]
	
	# save biom_table
	biom_name = tree_path.replace("tree.nwk", pf_name + "_table.txt")
	print(biom_name)
	biom_table.to_csv(biom_name, sep='\t')
