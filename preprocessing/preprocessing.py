import pandas as pd
import numpy as np
import glob
import os
import warnings

GFF_COLUMNS = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']
ATTRIBUTE_COLUMNS = ['ID', 'Parent', 'eC_number', 'Name', 'dbxref', 'gene', 
                     'inference', 'locus_tag', 'product', 'protein_id']
BED_COLUMNS = ['seqname', 'start', 'end', 'name']

# preprocessing combines multiple steps from previous pipelines 
# optional step 1: read the gff file 
def _read_annotation(path_to_annot):
    with open(path_to_annot) as f:
            lines = f.readlines()
    
    # check if FASTA sequence present; if so, get rid of it
    if '##FASTA' in lines:
        lines = lines[:lines.index('##FASTA')]
    
    # join the lines together into a string
    text_str = ''.join(lines)

    # turn each individual line into an element in a list 
    list_of_lines = text_str.split('\n')
    
    # remove comment lines; make tab-separated elements into new array
    array = []
    for line in list_of_lines:
        if (not line.startswith('#')) & (line != '') :
            array.append(line.split('\t'))

    # check that array has the correct number of columns 
    if array == [[]]: 
        warnings.warn(path_to_annot + 'does not contain genomic annotation information')
    elif np.array(array).shape[1] != 9: 
        raise ValueError(path_to_annot + 'does not have correct GFF dimensions')

    return pd.DataFrame(columns=GFF_COLUMNS, data=array)

# step 1: concatenate annotations together 
def concat_annotations(file_dict = None, glob_pattern=None): 
    """This function takes directory containing several 
    .gff/.csv files and concatenates them into a pd.DataFrame
    to make parsing easier for downstream steps.
    files_dir: directory containing .gff files 
    files_pattern: pattern to glob the .gff files 
    files_cols: ending column names desired
    to_remove: names of problematic files 
    """
    if (file_dict is None) & (glob_pattern is None):
        raise ValueError('One of file dictionary or path to files is required.')
    elif (file_dict is not None) & (glob_pattern is not None): 
        raise ValueError('Only one of file dictionary or path to files can be supplied.')
    
    dataframes = []

    if file_dict is not None:
        for key in file_dict.keys(): 
        # add the key to each dataframe as the filename column
            df = file_dict[key]
            df.insert(loc=0, column='filename', value=key)
            dataframes.append(df)

    if glob_pattern is not None:
        files = glob.glob(glob_pattern)
        if files == []:
            raise ValueError('Path contains no files')
        for f in files:
            if ('.gff' not in f) & ('.txt' not in f) & ('.tsv' not in f):
                raise ValueError('Path contains unexpected filetype')

            gff_read = _read_annotation(f)
            # extract filename and add it to the dataframe
            name = os.path.basename(f)
            end_index = name.index('.') # get rid of .gff 
            gff_read.insert(loc=0, column='filename', value=name[:end_index])
            dataframes.append(df)
    # concat everything in the list
    return pd.concat(dataframes).reset_index(drop=True)

# step 2: filter to coding sequences (or something else)
def filter_annotations(annotation_df, feature_value):
    """Filter annotations to a specific feature type (i.e. coding
    sequence, tRNAs, etc."""
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

def _split_attribute(df, cols_to_insert):
    """Helper function to split the attribute column in gff 
    files so that more meaningful information can be extracted"""
    # insert columns 
    for col_name in cols_to_insert:
        sub = col_name + '='
        df.insert(loc=0, column=col_name,
                  value=df['attribute'].apply(_sub_col,
                                              str_to_find=sub,
                                              sep=';'))
    # remove attribute
    df.drop(columns=['attribute'], inplace=True)
    return df

def extract_pfam(features_df):
    """Function extract Pfam features."""

    # split up the attribute column to search more easily for Pfams
    full_df = _split_attribute(features_df, ATTRIBUTE_COLUMNS)
    
    # find rows that have a Pfam value and insert column with their value 
    pfam = full_df.loc[full_df['inference'].str.contains('Pfam')]
    pfam.insert(loc=0, column='Pfam', 
                value=pfam['inference'].apply(_sub_col, str_to_find='Pfam-A:', 
                                              sep=':'))
    return pfam

# step 4: filter the features to certain value counts 
def filter_features(features_df, feature_col, value):
    counts = features_df[feature_col].value_counts()
    filt = counts[counts > value].index
    return features_df.loc[features_df[feature_col].isin(filt)]

# step 5: wrap to make a list of total features 
def wrapper_func(files_dir, files_pattern, pfam=True,
                 feature_value='CDS', filter_value=5):
    
    # get all annotations 
    annots = concat_annotations(files_dir, files_pattern)

    # filter annots to feature of interest
    feats = filter_annotations(annots, feature_value)

    # extract Pfams if flag is raised and filter Pfams
    if pfam: 
        feats = extract_pfam(feats)
        filtered = filter_features(feats, 'Pfam', filter_value)
    else: 
        filtered = filter_features(feats, 'gene', filter_value)
    
    return filtered
    
# step 6: make bed files and extract sequence
def _make_fasta_name(x, col_to_sort):
    return x[col_to_sort] + '_' + x['filename'] + '_' + str(x['rank'])

def extract_sequence(feats_df, fasta_dir, out_path, pfam=True):
    """Function to output bed files of all annotated pfams/genes"""
    # assign sort variable based on pfam flag
    if pfam: 
        sort='Pfam'
    else: 
        sort='gene'

    # apply groupby and apply functions based on sort variable    
    grouping = feats_df.groupby([sort, 'filename'])
    fasta_name = feats_df.apply(_make_fasta_name, args=(sort, ), axis=1)
    
    # insert rank and fasta_name columns
    feats_df.insert(loc=0, column='rank', 
                    value=grouping['start'].rank().astype('int'))
    feats_df.insert(loc=0, column='name', value=fasta_name)

    for f_id in feats_df[sort].unique():
        df = feats_df.loc[feats_df[sort]==f_id]
        unique_files = df['filename'].unique()

        if not os.path.exists(out_path + f_id):
            os.mkdir(out_path + f_id + '/')
        
        for f in unique_files:
            sub_df = df.loc[df['filename']==f]
            bed_name = out_path + f_id + '/' + f + '.bed'
            
            # ENH: use pybedtools instead of relying on CLI
            #bed_file = BedTool.from_dataframe(sub_df[BED_COLUMNS])
            #fasta_path = fasta_dir + f + '.fna'
            #subseq = bed_file.sequence(fi=fasta_path)
            
            # save the subsetted path 
            sub_df[BED_COLUMNS].to_csv(bed_name, index=False, header=False, sep='\t')

# in unix: concatenate fasta files and do msa
# def: convert sth to fasta 
# in unix: make all letters uppercase
# def: import msa into q2-phylogeny and run MAFFT 
# def: make feature table 
# in unix: convert txt to biom table