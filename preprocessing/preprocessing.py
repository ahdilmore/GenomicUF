import pandas as pd
import numpy as np
import glob
import os
import warnings
from typing import Tuple

GFF_COLUMNS = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']
ATTRIBUTE_COLUMNS = ['ID', 'Parent', 'eC_number', 'Name', 'dbxref', 'gene', 
                     'inference', 'locus_tag', 'product', 'protein_id']

# preprocessing combines multiple steps from previous pipelines 
# optional step 1: read the gff file 
def _read_annotation(path_to_annot) -> pd.DataFrame:
    with open(path_to_annot) as f:
            lines = f.readlines()
        
    # join the lines together into a string
    text_str = ''.join(lines)
    
    # check if FASTA sequence present; if so, get rid of it
    if '##FASTA' in text_str:
        text_str = text_str[:text_str.index('##FASTA')]

    # turn each individual line into an element in a list 
    list_of_lines = text_str.split('\n')
    
    # remove comment lines; make tab-separated elements into new array
    array = []
    for line in list_of_lines:
        if (not line.startswith('#')) & (line != '') :
            array.append(line.split('\t'))

    # check that array has the correct number of columns 
    if (array == [[]]) | (array==[]): 
        warnings.warn(path_to_annot + ' does not contain genomic annotation information')
        return pd.DataFrame([])
    return pd.DataFrame(columns=GFF_COLUMNS, data=array)

def process_data_dict(glob_pattern: str = None, 
                      data_dict : dict = None,
                      gff_ext : str = '.gff', 
                      fa_ext : str ='.fa') -> dict:
    # check that one or the other is present 
    if (data_dict is None) & (glob_pattern is None):
        raise ValueError('One of file dictionary or path to files is required.')
    elif (data_dict is not None) & (glob_pattern is not None): 
        raise ValueError('Only one of file dictionary or path to files can be supplied.')
    
    # check that the data dict is valid 
    if data_dict is not None:
        for key in data_dict.keys(): 
            # check that each input is a tuple of filepaths 
            if len(data_dict[key]) != 2:
                raise ValueError('More than two filepaths given in a data dictionary.')
            # check that each input is a valid path
            elif (gff_ext not in data_dict[key][0]):
                raise ValueError('Filepath given for annotations does not have ' + gff_ext + ' extension.')
            elif not (os.path.isfile(data_dict[key][0])):
                raise ValueError('The input ' + gff_ext + ' path does not exist.')
            elif (fa_ext not in data_dict[key][1]):
                raise ValueError('Filepath given for fasta argument does not have ' + fa_ext + ' extension.')
            elif not (os.path.isfile(data_dict[key][1])):
                raise ValueError('The input ' + fa_ext + ' path does not exist.')
        return data_dict

    elif glob_pattern is not None:
        directory_name = os.path.dirname(glob_pattern)
        sample_names = [os.path.basename(x) for x in glob.glob(glob_pattern)]
        data_dict = {}

        if len(sample_names) == 0: 
            raise ValueError('No files found after globbing!')
        
        for s in sample_names:
            # get fasta and gff
            gff_file = glob.glob(directory_name + '/' + s + '/*' + gff_ext)
            fasta_file = glob.glob(directory_name + '/' + s + '/*' + fa_ext)

            # check that these are unique 
            if len(gff_file) == 0: 
                raise ValueError('No ' + gff_ext + ' files found in a subdirectory.')
            elif len(fasta_file) == 0:
                raise ValueError('No ' + fa_ext + ' files found in a subdirectory.')
            elif len(gff_file) != 1:
                raise ValueError('More than one ' + gff_ext  + ' file found in a subdirectory.')
            elif len(fasta_file) != 1: 
                raise ValueError('More than one ' + fa_ext + ' file found in a subdirectory.')

            # add to data dictionary
            data_dict[s] = (gff_file[0], fasta_file[0])
        return data_dict

# step 1: concatenate annotations together 
def concat_annotations(data_dict : dict) -> pd.DataFrame: 
    '''Concatenates all annotations in the data dictionary together.'''
    dataframes = []
    for key in data_dict.keys(): 
        # read the gff file  
        df = _read_annotation(data_dict[key][0])
        # add the key to each dataframe as the filename column
        df.insert(loc=0, column='filename', value=key)
        dataframes.append(df)
    # concat everything in the list 
    return pd.concat(dataframes).reset_index(drop=True)

# step 2: filter to coding sequences (or something else)
def filter_annotations(annotation_df : pd.DataFrame, feature_value : str) -> pd.DataFrame:
    """Filter annotations to a specific feature type (i.e. coding
    sequence, tRNAs, etc."""
    if 'feature' not in annotation_df.columns:
        raise ValueError('The feature column is not found in the annotation dataframe.')
    elif feature_value not in annotation_df['feature'].unique():
        raise ValueError('The feature value provided is not present in the annotation dataframe.')
    return annotation_df.loc[annotation_df['feature']==feature_value]

# step 3 filter to features that have a Pfam annotation
def _sub_col(x, str_to_find, sep):
    """Helper function for extracting useful metrics from within a 
    string / column of pd.DataFrame"""
    if str_to_find in x:
        start_index = x.index(str_to_find)
        sub = x[start_index + len(str_to_find):]
        if sep in sub:
            end_index = sub.index(sep)
            return sub[:end_index]
        else:
            return sub
    else:
        return np.nan

def _split_attribute(df, cols_to_insert) -> pd.DataFrame:
    """Helper function to split the attribute column in gff 
    files so that more meaningful information can be extracted"""
    # insert columns 
    for col_name in cols_to_insert:
        sub = col_name + '='
        df[col_name] = df['attribute'].apply(_sub_col,
                                             str_to_find=sub,
                                             sep=';')
    # remove attribute
    return df.drop(columns=['attribute'])

def extract_pfam(features_df : pd.DataFrame) -> pd.DataFrame:
    """Function extract Pfam features."""
    if 'attribute' not in features_df.columns: 
        raise ValueError('Attribute column not present in dataframe.')

    # split up the attribute column to search more easily for Pfams
    full_df = _split_attribute(features_df, ATTRIBUTE_COLUMNS)

    # check that there are no NaNs in inference 
    if full_df.loc[full_df['inference'].isna()].shape[0] > 0:
        raise ValueError('Features dataframe has not been filtered. There are entries with invalid inferences.')
    
    # find rows that have a Pfam value and insert column with their value 
    pfam = full_df.loc[full_df['inference'].str.contains('Pfam')]
    pfam.insert(loc=0, column='Pfam', 
                value=pfam['inference'].apply(_sub_col, str_to_find='Pfam-A:', 
                                              sep=':'))
    return pfam

# step 4: filter the features to certain value counts 
def filter_features(features_df : pd.DataFrame, feature_col : str, min_value : int) -> pd.DataFrame:
    if (feature_col not in features_df.columns):
        raise ValueError('Column not found in dataframe provided.')
    counts = features_df[feature_col].value_counts()
    filt = counts[counts > min_value].index
    return features_df.loc[features_df[feature_col].isin(filt)]

# step 5: wrap to make a list of total features 
def wrapper_func(data_dict : dict = None, 
                 glob_pattern : str = None, 
                 gff_ext : str = '.gff', 
                 fa_ext : str = '.fa', 
                 pfam : bool = True, 
                 feature_value : str = 'CDS', 
                 filter_value : int =5) -> Tuple[dict, pd.DataFrame]:
    # check or make data dict 
    data_dict = process_data_dict(glob_pattern, data_dict, gff_ext, fa_ext)

    # get all annotations 
    annots = concat_annotations(data_dict)

    # filter annots to feature of interest
    feats = filter_annotations(annots, feature_value)

    # extract Pfams if flag is raised and filter Pfams
    if pfam: 
        pfam_feats = extract_pfam(feats)
        filtered = filter_features(pfam_feats, 'Pfam', filter_value)
    else: 
        filtered = filter_features(feats, 'gene', filter_value)
    
    return data_dict, filtered