import numpy as np
import pandas as pd
import sys
import pytest
sys.path.append('../')
from preprocessing.preprocessing import *

GFF_COLUMNS = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']
VALID_PATH = 'data/adaptation_AZ20/*S*'
VALID_DATA_DICT, VALID_OUT = wrapper_func(glob_pattern=VALID_PATH, fa_ext='.fna')
VALID_DF = concat_annotations(VALID_DATA_DICT)

def test_process_data_dict():
    # test neither data_dict nor glob_pattern given 
    with pytest.raises(ValueError, match='One of'):
        process_data_dict(fa_ext='.fna')
    # test both data_dict and glob_pattern given
    with pytest.raises(ValueError, match='Only one'):
        process_data_dict(glob_pattern=VALID_PATH, data_dict=VALID_DATA_DICT)
    
    # test invalid data_dicts 
    invalid_number = {'sample_1': ('/path/to/sample1.gff', '/path/to/sample1.fna', 'extra.fna'),
                      'sample_2': ('/path/to/sample2.gff', '/path/to/sample2.fna'), 
                      'sample_3': ('/path/to/sample3.gff', '/path/to/sample3.fna')}
    invalid_gff = {'sample_2': ('/path/to/sample2.gff', '/path/to/sample2.fna'), 
                   'sample_3': ('/path/to/sample3.gff', '/path/to/sample3.fna')}
    invalid_fa = {'sample_n': ('data/adaptation_AZ20/25N_Day7_10_S814/25N_Day7_10_S814_AZ20_mut.gff', 
                               'data/adaptation_AZ20/25N_Day7_10_S814/25N_Month3_10_S54_AZ20_mut.fna')}
    with pytest.raises(ValueError, match='More than two'):
        process_data_dict(data_dict=invalid_number)
    with pytest.raises(ValueError, match='.tsv extension'):
        process_data_dict(data_dict=invalid_gff, gff_ext='.tsv')
    with pytest.raises(ValueError, match='.gff path does not exist'):
        process_data_dict(data_dict=invalid_gff)
    with pytest.raises(ValueError, match='fasta argument'):
        process_data_dict(data_dict=invalid_fa)
    with pytest.raises(ValueError, match='.fna path does not exist'):
        process_data_dict(data_dict=invalid_fa, fa_ext='.fna')

    # test the invalid glob patterns
    no_fa = 'data/glob_paths/no_fa/*' 
    no_gff = 'data/glob_paths/no_gff/*'
    two_fa = 'data/glob_paths/extra_fa/*'
    two_gff = 'data/glob_paths/extra_gff/*'
    no_files = 'data/glob_paths/wrong_dir/*'
    with pytest.raises(ValueError, match='No files found'):
        process_data_dict(glob_pattern=no_files)
    with pytest.raises(ValueError, match='No .gff files'):
        process_data_dict(glob_pattern=no_gff, fa_ext='.fna')
    with pytest.raises(ValueError, match='No .fna files'):
        process_data_dict(glob_pattern=no_fa, fa_ext='.fna')
    with pytest.raises(ValueError, match='More than one .gff'):
        process_data_dict(glob_pattern=two_gff, fa_ext='.fna')
    with pytest.raises(ValueError, match='More than one .fna'):
        process_data_dict(glob_pattern=two_fa, fa_ext='.fna')
    
    valid_pattern = 'data/glob_paths/valid/*'
    valid_out = {'sample_2': ('data/glob_paths/valid/sample_2/25N_Day7_10_S814_AZ20_mut.gff',
                              'data/glob_paths/valid/sample_2/25N_Day7_10_S814_AZ20_mut.fna'),
                 'sample_1': ('data/glob_paths/valid/sample_1/46R_Day7_13_S292_AZ51_mut.gff',
                              'data/glob_paths/valid/sample_1/46R_Day7_13_S292_AZ51_mut.fna')}
    assert process_data_dict(glob_pattern=valid_pattern, fa_ext='.fna') == valid_out

def test_concat_annotations():
    # check that the warning shows up 
    az51_path = 'data/adaptation_AZ51/*S*'
    az51_data_dict = process_data_dict(glob_pattern=az51_path, fa_ext='.fna')
    with pytest.warns(UserWarning, match='does not contain'):
        concat_annotations(data_dict=az51_data_dict)
    az51_df = concat_annotations(data_dict=az51_data_dict)
    assert az51_df.shape == (35432, 10)
 
def test_filter_annotations():
    invalid_df = VALID_DF.rename(columns={'feature': 'feature_id'})
    # check that feature column exists in dataframe
    with pytest.raises(ValueError, match='column is not found'):
        filter_annotations(annotation_df=invalid_df, feature_value='CDS')
    # check that value exists in col 
    with pytest.raises(ValueError, match='value provided'):
        filter_annotations(annotation_df=VALID_DF, feature_value='mRNA')
    cds_only = filter_annotations(annotation_df=VALID_DF, feature_value='CDS')
    assert cds_only.shape == (28907, 10)

# test pfam filter 
def test_extract_pfam(): 
    # check that the inference column is present in the dataframe 
    invalid_df = VALID_DF.drop(columns=['attribute'])
    with pytest.raises(ValueError, match='Attribute column'):
        extract_pfam(invalid_df)
    # when there are invalid inferences, value error should be raised
    with pytest.raises(ValueError, match='invalid inferences'):
        extract_pfam(VALID_DF)
    valid_cds = filter_annotations(annotation_df=VALID_DF, feature_value='CDS')
    pfam_valid = extract_pfam(valid_cds)
    assert pfam_valid.shape == (1872, 20)

def test_filter_features(): 
    # check that col exists in dataframe 
    with pytest.raises(ValueError, match='Column not found'):
        filter_features(features_df=VALID_DF, feature_col='feature_id', min_value=5)
    valid_cds = filter_annotations(annotation_df=VALID_DF, feature_value='CDS')
    pfam_valid = extract_pfam(valid_cds)
    feat_filt = filter_features(features_df=pfam_valid, feature_col='Pfam', min_value=10)
    assert feat_filt.shape == (199, 20)