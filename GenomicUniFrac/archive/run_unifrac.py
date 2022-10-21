import pandas as pd
import numpy as np
import glob
import skbio
import unifrac
import biom
from skbio import TreeNode

fa_dir='genes_out_AZ51_filt/'
md_df = pd.read_csv("AZ51_full_metadata.csv", index_col="sample_name")
tree_files=glob.glob(fa_dir + '*tree/tree.nwk')

col_names = ['n_samples', 'n_groups', 'pseudo_f', 'p_val']
base_df = pd.DataFrame(columns=col_names)

for tree_path in tree_files:
	biom_path = tree_path.replace("tree.nwk", "table.biom")
	table = biom.load_table(biom_path)
	tree = skbio.io.read(tree_path, format='newick', into=TreeNode)	

	# figure out the name 
	pf_index = tree_path.find("AZ51_filt") + 10
	tree_index = tree_path.find("_tree")
	pf_name = tree_path[pf_index:tree_index]	

	# subset the metadata
	md_subset = md_df.loc[table.ids()] 
	
	if(len(md_subset['time_point'].unique()) == 1): 
		print(pf_name + ' only one time point')
		data_pd = {col_names[0]:table.shape[1], col_names[1]:1,
			   col_names[2]:np.nan, col_names[3]:np.nan}
	elif(tree.count() == tree.count(tips=True)+1):
		print(pf_name + ' no internal nodes')
		data_pd = {col_names[0]:table.shape[1],
			   col_names[1]:len(md_subset['time_point'].unique()),
			   col_names[2]:np.nan, col_names[3]:np.nan}
	else:
		uwuf_dm = unifrac.unweighted(table=biom_path, phylogeny=tree_path)
		uwuf_permanova = skbio.stats.distance.permanova(uwuf_dm, md_subset,
							        "time_point",
								permutations=10)
		data_pd = {col_names[0]:uwuf_permanova[2], col_names[1]:uwuf_permanova[3],
                           col_names[2]:uwuf_permanova[4], col_names[3]:uwuf_permanova[5]}
	base_df = pd.concat([base_df, pd.DataFrame(data=data_pd, index=[pf_name])])
	base_df.to_csv('gene_results_AZ51.csv')	
