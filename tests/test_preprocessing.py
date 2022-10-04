import numpy as np
import pandas as pd
import pytest
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

