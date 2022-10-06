import numpy as np
import pandas as pd
import pytest
# gonna have to figure out why this is giving me an error
from preprocessing.preprocessing import *

GFF_COLUMNS = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']

# tests for concatenating dataframes
def test_invalid_paths(): 
    # no files present 
    invalid_path = 'data/hazels_files/'
    with pytest.raises(ValueError, match='no files'):
        concat_annotations(glob_pattern=invalid_path)

    # files present, but has some besides gffs 
    no_gffs = 'data/adaptation_AZ51/46R_Day7_13_S292/'
    with pytest.raises(ValueError, match='unexpected filetype'):
        concat_annotations(glob_pattern=no_gffs)

def test_valid_files():
    valid_path = 'data/adaptation_AZ20/*S*/*.gff'
    # show that no errors are raised
    valid_out = concat_annotations(glob_pattern=valid_path)
    assert valid_out.shape == (58906, 10)

def test_empty_file_warning():
    # checks that user is warned about files without genomic information
    contains_empty_file = 'data/adaptation_AZ51/*S*/*.gff'
    with pytest.warns(UserWarning, match='does not contain'):
        concat_annotations(glob_pattern=contains_empty_file)

def test_gff_columns_present(): 
    # data setup 
    data_1 = [[1, 2, 3, 4, 5, 6, 7, 8, 9],
              [10, 11, 12, 13, 14, 15, 16, 17, 18]]
    data_2 = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]]
    data_3 = [[9, 8, 7, 6, 5, 4, 3, 2, 1]]
    incorrect_columns = {
        'file_1' : pd.DataFrame(columns=GFF_COLUMNS, data=data_1), 
        'file_2' : pd.DataFrame(data=data_2)
    }
    incorrect_col_names = {
        'file_1' : pd.DataFrame(columns=GFF_COLUMNS, data=data_1), 
        'file_3' : pd.DataFrame(data=data_3)
    }
    correct_names = {
        'file_1' : pd.DataFrame(columns=GFF_COLUMNS, data=data_1), 
        'file_3' : pd.DataFrame(columns=GFF_COLUMNS, data=data_3)
    }

    # input dataframes have incorrect column dimensions
    with pytest.raises(ValueError, match='dimensions'):
        concat_annotations(file_dict=incorrect_columns)
    
    with pytest.raises(ValueError, match='columns'):
        concat_annotations(file_dict=incorrect_col_names)
    
    # check that correct_names gives no error
    concat_correct_format = concat_annotations(file_dict=correct_names)
    assert sum(concat_correct_format.columns == ['filename'] + GFF_COLUMNS) == 10
    assert concat_correct_format.shape == (3, 10)

# test filter_annotations 
def test_filter_annotations():
    valid_path = 'data/adaptation_AZ20/*S*/*.gff'
    valid_df = concat_annotations(glob_pattern=valid_path)
    invalid_df = valid_df.rename(columns={'feature': 'feature_id'})
    # check that feature column exists in dataframe
    with pytest.raises(ValueError, match='column is not found'):
        filter_annotations(annotation_df=invalid_df, feature_value='CDS')
    # check that value exists in col 
    with pytest.raises(ValueError, match='value provided'):
        filter_annotations(annotation_df=valid_df, feature_value='mRNA')
    cds_only = filter_annotations(annotation_df=valid_df, feature_value='CDS')
    assert cds_only.shape == (28907, 10)

# test pfam filter 
def test_pfam_filter(): 
    # check that the inference column is present in the dataframe 
    valid_path = 'data/adaptation_AZ20/*S*/*.gff'
    valid_df = concat_annotations(glob_pattern=valid_path)
    invalid_df = valid_df.drop(columns=['attribute'])
    with pytest.raises(ValueError, match='Attribute column'):
        extract_pfam(invalid_df)
    # check that the output is as expected in valid dataframe 
    # having some issues with this check but i'm hungry
    pfams_valid = extract_pfam(valid_df)
    #assert pfams_valid.shape

def test_filter_features(): 
    valid_path = 'data/adaptation_AZ20/*S*/*.gff'
    valid_df = concat_annotations(glob_pattern=valid_path)
    # check that col exists in dataframe 
    with pytest.raises(ValueError, match='Column not found'):
        filter_features(features_df=valid_df, feature_col='feature_id', value=5)
    # check that min value is an integer
    with pytest.raises(ValueError, match='integer'):
        filter_features(features_df=valid_df, feature_col='feature', value='CDS')
    feat_filt = filter_features(feature_df=valid_df, feature_col='Pfam', value=10)
    # assert correct dims 