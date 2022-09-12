import pandas as pd
import itertools
import unifrac
import skbio

# inputs
df = pd.read_csv("pfam_results_AZ51.csv", index_col='Unnamed: 0')
md = pd.read_csv("AZ51_metadata.csv", index_col="sample_name")
input_dir = 'out_AZ51/msa_files/'
col_names = ['PF_1', 'PF_2', 'n_samples', 'n_groups', 'pseudo_f', 'p_val']
base_df = pd.DataFrame(columns=col_names)

# identify top 100 pfams
top100_df = df.loc[df["n_samples"] > 50].sort_values(by="pseudo_f", ascending=False)[:100]

# make list of combinations 
for combo in itertools.combinations(top100_df.index, 2):
	dir1 = input_dir + combo[0] + '_tree/'
	dir2 = input_dir + combo[1] + '_tree/'
	meta_uni = unifrac.meta(tables=[dir1 + 'table.biom', dir2 + 'table.biom'], 
				phylogenies=[dir1 + 'tree.nwk', dir2 + 'tree.nwk'],
				method='unweighted')
	meta_permanova = skbio.stats.distance.permanova(meta_uni, md, "time_point",
							permutations=1)
	data_pd = {col_names[0]:combo[0], col_names[1]:combo[1],
		   col_names[2]:meta_permanova[2], col_names[3]:meta_permanova[3],
		   col_names[4]:meta_permanova[4], col_names[5]:meta_permanova[5]}
	combo_name = ','.join(list(combo))
	base_df = pd.concat([base_df, pd.DataFrame(data=data_pd, index=[combo_name])])
	base_df.to_csv('meta_top_results_AZ51.csv')
	print(combo_name)
