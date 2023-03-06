import pandas as pd
import numpy as np 
import glob
import skbio
import os
import unifrac
import itertools
import biom
    
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

def single_gene(unifracs_to_run : list,
                metadata : pd.DataFrame, 
                sep_column : str, 
                table : biom.Table = None, 
                tree_dir : str = None, 
                table_and_tree_dir : str = None):
    # Checks that either table AND tree dir OR table_and_tree_dir are passed
    if ((tree_dir is None) | (table is None)) & (table_and_tree_dir is None):
        raise ValueError('Either a directory containing .biom tables and .nwk trees or a .biom table and directory with .nwk trees must be passed.')
    
    # if table and tree dir provided, tree glob added 
    if (table_and_tree_dir is not None): 
        tree_dir = table_and_tree_dir
    # Checks to make sure all unifracs_to_run are viable inputs
    for element in unifracs_to_run:
        if element not in UNIFRACS:
            raise ValueError(element + " is not a valid unifrac input")
    
    names = []
    test_stats = []
    p_values = []
    unifrac_type = []
    
    for unifrac_method in unifracs_to_run:
        for path in glob.glob(tree_dir + "*.nwk")[:2]:
            if len(os.path.basename(path).split('.')) < 2:
                raise ValueError(os.path.basename(path) + " is not a valid tree file name")
            if table is None:
                table = biom.load_table(path.replace('tree.nwk', 'table.biom'))
                # filter table based on metadata
                table.filter(ids_to_keep=metadata.index.intersection(table.ids())) 
            tree = skbio.io.read(path, format="newick", into=skbio.TreeNode)

            names.append(os.path.basename(path).split('.')[-2])
            perm_out = run_unifracs(table, tree, metadata, sep_column, UNIFRACS[unifrac_method])
            test_stats.append(perm_out['test statistic'])
            p_values.append(perm_out['p-value'])
            unifrac_type.append(unifrac_method)
            
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
    
    names = []
    test_stats = []
    p_values = []
    unifrac_type = []
    
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