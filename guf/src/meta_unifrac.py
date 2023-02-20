import pandas as pd
import numpy as np 
import glob
import skbio
import os
import unifrac
import biom
import itertools

table = biom.load_table('Data/thdmi_feature_table.biom')
metadata = pd.read_csv('Data/consolidated_metadata.tsv', sep='\t',index_col=0)

    
# add any additional viable unifrac methods here
UNIFRACS = {
    'unweighted': unifrac.unweighted,
    'unweighted fp32': unifrac.unweighted_fp32,
    'meta': unifrac.meta
}

def run_unifracs(table, tree, metadata, column, unifracToRun, method='None'):
    if method == 'None':
        dm = unifracToRun(table,tree)
    else:
        dm = unifracToRun(table,tree,method=method)
    return skbio.stats.distance.permanova(dm, metadata, column)

# can add a threads parameter later if we want to

def single_gene(tree_dir, unifracs_to_run, table, metadata, sep_column):
    
    valid_unifracs = []
    
    # Checks to make sure all unifracs_to_run are viable inputs
    for element in unifracs_to_run:
        if element not in UNIFRACS.keys():
            warnings.warn(element + " is not a valid unifrac input")
        else:
            valid_unifracs.append(element)
    
    print("Unifracs are Valid")
    
    names = []
    test_stats = []
    p_values = []
    unifrac_type = []
    
    i=0
    for unifrac_method in unifracs_to_run:
        for path in glob.glob(tree_dir + "*.nwk")[:2]:
            if len(os.path.basename(path).split('.')) < 2:
                warnings.warn(os.path.basename(path) + " is not a valid tree file name")
            names.append(os.path.basename(path).split('.')[-2])
            
            tree = skbio.io.read(path, format="newick", into=skbio.TreeNode)
            perm_out = run_unifracs(table, tree, metadata, sep_column, UNIFRACS[unifrac_method])
            test_stats.append(perm_out['test statistic'])
            p_values.append(perm_out['p-value'])
            unifrac_type.append(unifrac_method)
            
            # This is just here for testing and such
            print(str(i + 1) + ".\t" + names[i] + " \t| " + str(test_stats[i]) + " \t| " + str(p_values[i]) + " \t| " + unifrac_type[i])
            
            i = i + 1
            
    df = pd.DataFrame(data = {'test statistic': test_stats,
                                     'p value': p_values,
                                     'unifrac type': unifrac_type},
                     index = names)
    
    return df

def multi_gene(tree_dir, tables, methods, metadata, sep_column):
    
    valid_unifracs = []
    
    # Checks to make sure all unifracs_to_run are viable inputs
    for element in methods:
        if element not in UNIFRACS.keys():
            warnings.warn(element + " is not a valid unifrac input")
        else:
            valid_unifracs.append(element)
    
    print("Unifracs are Valid")
    
    names = []
    test_stats = []
    p_values = []
    unifrac_type = []
    
    i=0
    
    combo_length = len(tables)
    
    for method in methods:
        for path_combo in list(itertools.combinations(glob.glob(tree_dir+'*.nwk'), combo_length))[:5]:
            trees = [skbio.io.read(path, format='newick', into=skbio.TreeNode) for path in path_combo]
            tree_names = [os.path.basename(path).split('.')[-2] for path in path_combo]
            names.append(tree_names)
            perm_out = run_unifracs(tables, trees, metadata, sep_column, UNIFRACS['meta'], method=method)
            test_stats.append(perm_out['test statistic'])
            p_values.append(perm_out['p-value'])
            unifrac_type.append('meta ' + method)
            print(str(i + 1) + ".\t" + str(names[i]) + " \t| " + str(test_stats[i]) + " \t| " + str(p_values[i]) + " \t| " + unifrac_type[i])
            
            i = i + 1
            
        names = [str(name_combo) for name_combo in names]
            
    df = pd.DataFrame(data = {'test statistic': test_stats,
                                     'p value': p_values,
                                     'unifrac type': unifrac_type},
                     index = names)
    
    return df



def wrangled(df, unifracs_to_run):
    # below is just data wrangling and reorganizing
    # should I just make this into a separate function
    
    df.rename_axis(index = 'names', inplace=True)
    df = df.reset_index()
    df = df.set_index(["names","unifrac type"])
    df = df.unstack(level="unifrac type")
    col_array = []
    for uf in unifracs_to_run:
        col_array.append(('test statistic', uf))
        col_array.append(('p value', uf))
    df = df[col_array]
    
    df = df.swaplevel(axis=1)
    
    return df
            
            
single_df = single_genes("PerGeneTrees/", ["unweighted","unweighted fp32"], table, metadata, "thdmi_cohort")
wrangled_single_df = wrangled(single_df, ["unweighted","unweighted fp32"])

meta_df = multi_gene("PerGeneTrees/", [table,table], ["unweighted"], metadata, "thdmi_cohort")
wrangled_meta_df = wrangled(meta_df, ["meta unweighted"])

# UNIFRAC_COLUMNS = ['n_samples', 'n_groups', 'pseudo_f', 'p_val']

# def _get_unifrac_output(tree_path, table, metadata, 
#                         col_of_interest, permutations):
#     tree = skbio.io.read(tree_path, format='newick', into=skbio.TreeNode)
#     feat_id = os.path.basename(os.path.dirname(tree_path))
#     if len(metadata[col_of_interest].unique()) == 1:
#         # warn the user that only one group has this given gene
#         data_pd = {UNIFRAC_COLUMNS[0]:table.shape[1], 
#                    UNIFRAC_COLUMNS[1]:1,
#                    UNIFRAC_COLUMNS[2]:np.nan,
#                    UNIFRAC_COLUMNS[3]:np.nan}
#     elif tree.count() == (tree.count(tips=True) + 1): 
#         # warn the user that this tree has no internal nodes
#         data_pd = {UNIFRAC_COLUMNS[0]:table.shape[1],
#                    UNIFRAC_COLUMNS[1]:len(metadata[col_of_interest].unique()),
#                    UNIFRAC_COLUMNS[2]:np.nan, 
#                    UNIFRAC_COLUMNS[3]:np.nan}
#     else:
#         # todo: add other unifrac types
#         dm = unifrac.unweighted(table = tree_path.replace("tree.nwk", "table.biom"),
#                                 phylogeny = tree_path)
#         permanova = skbio.stats.distance.permanova(dm, metadata, col_of_interest, 
#                                                    permutations)
#         data_pd = {UNIFRAC_COLUMNS[0]:permanova[2], 
#                    UNIFRAC_COLUMNS[1]:permanova[3],
#                    UNIFRAC_COLUMNS[2]:permanova[4], 
#                    UNIFRAC_COLUMNS[3]:permanova[5]}
#     return pd.DataFrame(data=data_pd, index=[feat_id])
        

# def single_genes(preprocessed_dir : str,
#                  col_of_interest : str, 
#                  tree_ext : str = 'TreeData/',
#                  metadata : pd.DataFrame = None, 
#                  path_to_metadata : str = None, 
#                  permutations : int = 999):
#     paths_to_trees = glob.glob(preprocessed_dir + tree_ext + '*/tree.nwk')

#     # check that one of metadata_ext or metadata file is present 
#     if (path_to_metadata is None) & (metadata is None):
#         raise ValueError('One of metadata or path to metadata is required.')
#     elif (path_to_metadata is not None) & (metadata is not None): 
#         raise ValueError('Only one of metadata or path to metadata can be supplied.')
    
#     # if metadata_ext given, load the metadata
#     if path_to_metadata is not None: 
#         metadata = pd.read_csv(path_to_metadata, sep='\t', index_col='sample_name')
#     # check that the sample names in the metadata match the sample names in the output dir 
    
#     list_of_dfs = []
#     for tree_path in paths_to_trees: 
#         # load tree and table info, as well as feature ID
#         table = biom.load_table(tree_path.replace("tree.nwk", "table.biom"))
        
#         # subset the metadata
#         md_subset = metadata.loc[table.ids()]

#         # go through each tree and append unifrac results 
#         list_of_dfs.append(_get_unifrac_output(tree_path, 
#                                                table, 
#                                                md_subset, 
#                                                col_of_interest, 
#                                                permutations))
#     return pd.concat(list_of_dfs)
