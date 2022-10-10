import pandas as pd
import os 
from pybedtools import BedTool
BED_COLUMNS = ['seqname', 'start', 'end', 'name']

def _make_fasta_name(x, col_to_sort):
    return x[col_to_sort] + '_' + x['filename'] + '_' + str(x['rank'])

def extract_sequence(data_dict, feats_df, col_to_sort, out_path):
    """Function to output bed files of all annotated pfams/genes"""
    # ensure that feats_df.start and feats_df.end are integers
    feats_df.start = feats_df.start.astype('int')
    feats_df.end = feats_df.end.astype('int')
    # get rank for tree construction 
    feats_grouped = feats_df.groupby([col_to_sort, 'filename'])
    feats_df.insert(loc=0, column='rank', 
                    value=feats_grouped['start'].rank().astype('int'))
    # make name of fasta file using ranks calculated above 
    fasta_names = feats_df.apply(_make_fasta_name, args=(col_to_sort, ), axis=1)
    feats_df.insert(loc=0, column='name', value=fasta_names)

    wd = os.getcwd()
    for f_id in feats_df[col_to_sort].unique():
        f_df = feats_df.loc[feats_df[col_to_sort]==f_id]
        unique_files = feats_df['filename'].unique()

        if not os.path.exists(out_path + f_id):
            os.mkdir(out_path + f_id + '/')
        
        for f in unique_files:
            sub_df = f_df.loc[f_df['filename']==f]
            bed_file = BedTool.from_dataframe(sub_df[BED_COLUMNS])
            fasta_path = data_dict[f][1]
            bed_file.sequence(fi=fasta_path, name=True, fo=out_path+f_id+'/'+f+'.fa')
        
        os.chdir(out_path + f_id)
        os.system("cat *.fa > merged.fa")
        os.chdir(wd)