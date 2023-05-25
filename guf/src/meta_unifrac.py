import pandas as pd
import numpy as np 
import glob
import skbio
import bp
import os
import unifrac
import itertools
import biom
    
# add any additional viable unifrac methods here
UNIFRACS = {
    'unweighted': unifrac.unweighted,
    'unweighted_fp32': unifrac.unweighted_fp32,
    'weighted': unifrac.weighted_unnormalized, 
    'weighted_fp32': unifrac.weighted_unnormalized_fp32,
    'meta': unifrac.meta
}

def _verify_unifracs(unifracs_to_run): 
    for element in unifracs_to_run:
        if element not in UNIFRACS:
            raise ValueError(element + " is not a valid unifrac input")

def run_unifracs(table, tree, metadata, column, unifracToRun, method=None, permutations=1):
    if method is None:
        dm = unifracToRun(table,tree)
        intersect = metadata.index.intersection(table.ids())
    else:
        dm = unifracToRun(table,tree,method=method)
        t = table[0].ids()
        for i in range(1, len(table)):
            t = set(t) | set(table[i].ids())
        intersect = metadata.index.intersection(list(t))
    dm_filt = dm.filter(ids=intersect)
    if len(metadata.loc[intersect][column].unique()) == 1: 
        return None
    return skbio.stats.distance.permanova(dm_filt, metadata, column, permutations)

def single_gene(unifracs_to_run : list,
                sample_metadata : pd.DataFrame, 
                sep_column : str,
                feature_metadata : pd.DataFrame = None,
                table : biom.Table = None, 
                tree_dir : str = None, 
                table_and_tree_dir : str = None):
    # Checks that either table AND tree dir OR table_and_tree_dir are passed
    if ((tree_dir is None) | (table is None)) & (table_and_tree_dir is None):
        raise ValueError('Either a directory containing .biom tables and .nwk trees or a .biom table and directory with .nwk trees must be passed.')
    
    # if table and tree dir provided, tree glob added 
    if (table_and_tree_dir is not None): 
        tree_dir = table_and_tree_dir
        table_dict = {}
    tree_dict = {}
    _verify_unifracs(unifracs_to_run)
    
    for path in glob.glob(tree_dir + "*.nwk"):
        tree_dict[path] = bp.parse_newick(open(path).read())
        if table_and_tree_dir is not None: 
            table_dict[path] = biom.load_table(path.replace('tree.nwk', 'table.biom'))

    all_results = []
    for unifrac_method in unifracs_to_run:
        results = pd.DataFrame(columns=['PERMANOVA_PseudoF', 'p_value'])
        for path in glob.glob(tree_dir + "*.nwk"):
            if len(os.path.basename(path).split('.')) < 2:
                raise ValueError(os.path.basename(path) + " is not a valid tree file name")
            if table_and_tree_dir is not None:
                table = table_dict[path]
                name = os.path.basename(os.path.dirname(path))
            else: 
                name = os.path.basename(path).split('.')[0]
            tree = tree_dict[path]
            perm_out = run_unifracs(table, tree, sample_metadata, sep_column, UNIFRACS[unifrac_method])
            if perm_out is not None:
                results = pd.concat([results, pd.DataFrame(data = {'PERMANOVA_PseudoF': perm_out['test statistic'], 
                                                                   'p_value': perm_out['p-value']}, index=[name])])
            else: 
                results = pd.concat([results, pd.DataFrame(data = {'PERMANOVA_PseudoF': np.nan, 
                                                                   'p_value': np.nan}, index=[name])])
        results['unifrac_type'] = unifrac_method
        if feature_metadata is not None: 
            results = results.merge(feature_metadata, right_index=True, left_index = True)
        all_results.append(results)
    return pd.concat(all_results).reset_index()

def multi_gene(unifracs_to_run : list, tree_dir : str, sample_metadata, sep_column, 
               num_tables : int, table : biom.Table = None, max_count : int = None, 
               subset_to_run : list = None):
        
    # Checks to make sure all unifracs_to_run are viable inputs
    _verify_unifracs(unifracs_to_run)
    
    # Checks that either table AND tree dir OR table_and_tree_dir are passed
    if table is not None: 
        tables = []
        for i in range(num_tables): 
            tables.append(table)
    else: 
        table_dict = {}
    if subset_to_run is None: 
        combos = list(itertools.combinations(glob.glob(tree_dir+'*.nwk'), num_tables))
    elif len(subset_to_run) <= num_tables:
        raise ValueError("Number of tables to compare is larger than number of tables provided")
    else: 
        combos = list(itertools.combinations(subset_to_run, num_tables))
    all_results = []
    
    counter = 0

    tree_dict = {}
    for path in glob.glob(tree_dir+'*.nwk'): 
        tree_name = os.path.basename(os.path.dirname(path))
        tree_dict[tree_name] = bp.parse_newick(open(path).read())
        if table is None: 
            table_dict[tree_name] = biom.load_table(path.replace('tree.nwk', 'table.biom'))

    for method in unifracs_to_run:
        results = pd.DataFrame(columns=['PERMANOVA_PseudoF', 'p_value'])
        for path_combo in combos:
            if subset_to_run is None: 
                path_combo = [os.path.basename(os.path.dirname(i)) for i in path_combo]
            trees = [tree_dict[path] for path in path_combo]
            if table is None: 
                tables = [table_dict[p] for p in path_combo]
            names = '&'.join(path_combo)
            perm_out = run_unifracs(tables, trees, sample_metadata, sep_column, UNIFRACS['meta'], method=method)
            if perm_out is not None: 
                results = pd.concat([results, pd.DataFrame(data = {'PERMANOVA_PseudoF': perm_out['test statistic'], 
                                                                   'p_value': perm_out['p-value']}, index=[names])])
            else: 
                results = pd.concat([results, pd.DataFrame(data = {'PERMANOVA_PseudoF': np.nan, 
                                                                   'p_value': np.nan}, index=[names])])
            counter += 1
            if counter == max_count: 
                break 
        results['unifrac_type'] = 'meta' + method
        all_results.append(results)

    return pd.concat(all_results).reset_index()



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
