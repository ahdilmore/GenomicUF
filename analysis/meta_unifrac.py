import pandas as pd
import numpy as np 
import glob
import skbio
import os
import unifrac
import biom

UNIFRAC_COLUMNS = ['n_samples', 'n_groups', 'pseudo_f', 'p_val']

def _get_unifrac_output(tree_path, table, metadata, 
                        col_of_interest, permutations):
    tree = skbio.io.read(tree_path, format='newick', into=skbio.TreeNode)
    feat_id = os.path.basename(os.path.dirname(tree_path))
    if len(metadata[col_of_interest].unique()) == 1:
        # warn the user that only one group has this given gene
        data_pd = {UNIFRAC_COLUMNS[0]:table.shape[1], 
                   UNIFRAC_COLUMNS[1]:1,
                   UNIFRAC_COLUMNS[2]:np.nan,
                   UNIFRAC_COLUMNS[3]:np.nan}
    elif tree.count() == (tree.count(tips=True) + 1): 
        # warn the user that this tree has no internal nodes
        data_pd = {UNIFRAC_COLUMNS[0]:table.shape[1],
                   UNIFRAC_COLUMNS[1]:len(metadata[col_of_interest].unique()),
                   UNIFRAC_COLUMNS[2]:np.nan, 
                   UNIFRAC_COLUMNS[3]:np.nan}
    else:
        # todo: add other unifrac types
        dm = unifrac.unweighted(table = tree_path.replace("tree.nwk", "table.biom"),
                                phylogeny = tree_path)
        permanova = skbio.stats.distance.permanova(dm, metadata, col_of_interest, 
                                                   permutations)
        data_pd = {UNIFRAC_COLUMNS[0]:permanova[2], 
                   UNIFRAC_COLUMNS[1]:permanova[3],
                   UNIFRAC_COLUMNS[2]:permanova[4], 
                   UNIFRAC_COLUMNS[3]:permanova[5]}
    return pd.DataFrame(data=data_pd, index=[feat_id])
        

def single_genes(preprocessed_dir : str,
                 col_of_interest : str, 
                 tree_ext : str = 'TreeData/',
                 metadata : pd.DataFrame = None, 
                 path_to_metadata : str = None, 
                 permutations : int = 999):
    paths_to_trees = glob.glob(preprocessed_dir + tree_ext + '*/tree.nwk')

    # check that one of metadata_ext or metadata file is present 
    if (path_to_metadata is None) & (metadata is None):
        raise ValueError('One of metadata or path to metadata is required.')
    elif (path_to_metadata is not None) & (metadata is not None): 
        raise ValueError('Only one of metadata or path to metadata can be supplied.')
    
    # if metadata_ext given, load the metadata
    if path_to_metadata is not None: 
        metadata = pd.read_csv(path_to_metadata, sep='\t', index_col='sample_name')
    
    list_of_dfs = []
    for tree_path in paths_to_trees: 
        # load tree and table info, as well as feature ID
        table = biom.load_table(tree_path.replace("tree.nwk", "table.biom"))
        
        # subset the metadata
        md_subset = metadata.loc[table.ids()]

        # go through each tree and append unifrac results 
        list_of_dfs.append(_get_unifrac_output(tree_path, 
                                               table, 
                                               md_subset, 
                                               col_of_interest, 
                                               permutations))
    return pd.concat(list_of_dfs)