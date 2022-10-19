import pytest
import sys
sys.path.append('../')
from analysis.meta_unifrac import *

def test_single_genes():
    path_to_metadata_az51 = '../tests/data/adaptation_AZ20/metadata.csv'
    path_to_metadata_az20 = '../tests/data/adaptation_AZ51/metadata.csv'
    
    # make as a dataframe
    time_points = ['Day7', 'Month1', 'Month3', 'Month6']
    lanes = ['S292', 'S1816', 'S200', 'S678']
    metadata_df = pd.DataFrame(data={'time_point': time_points, 'sequencing_lane': lanes})
    metadata_df['mouse_name'] = '46R'
    metadata_df['isolate_number'] = 13
    metadata_df['bacteria'] = 'AZ51'
    metadata_df['sample_name'] = (metadata_df['mouse_name'] + '_' + 
                                  metadata_df['time_point'] + '_' + 
                                  metadata_df['isolate_number'].astype('str') + '_' + 
                                  metadata_df['sequencing_lane'])
    metadata_df = metadata_df.set_index('sample_name')

    with pytest.raises(ValueError, match='One of'):
        single_genes(preprocessed_dir='out_AZ20/', col_of_interest='time_point')
    with pytest.raises(ValueError, match='Only one'):
        single_genes(preprocessed_dir='out_AZ20/', col_of_interest='time_point', 
                     path_to_metadata=path_to_metadata_az20, 
                     metadata=metadata_df)
    # todo: test where the sample names in the metadata don't match the sample names in the tree
    with pytest.raises(ValueError, match='sample names'):
        single_genes(preprocessed_dir='out_AZ20/', col_of_interest='time_point', 
                     metadata=metadata_df)