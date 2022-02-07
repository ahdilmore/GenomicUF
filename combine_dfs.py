import pandas as pd
import numpy as np
import glob

start_dir = '../../callaban/Projects/Adaptation/prokka_AZ20/'
files = glob.glob(start_dir + '*L001/*_mut.gff')
final_cols = ['filename', 'seqname', 'source', 'feature', 'start',
              'end', 'score', 'strand', 'frame', 'attribute']
dataframes = []
for f in files:
	# open files
        with open(f) as g:
                lines = g.readlines()

        # find the index of fasta and cut off lines there
        lines = lines[:lines.index('##FASTA\n')]

        # join the string together
        big_str = ''.join(lines)
        list_of_lines = big_str.split('\n')

        # split line by tabs; remove comment lines
        array = []
        for line in list_of_lines:
                if line.startswith('#') == False:
                        array.append(line.split('\t'))

        # put into dataframe    
        df = pd.DataFrame(columns = ['seqname', 'source', 'feature',
                                     'start', 'end', 'score',
                                     'strand', 'frame', 'attribute'],
                          data = array)
        df = df.dropna()
        name = f[len(start_dir):]
        start_index = name.index('/') + 1
        end_index = name.index('.')
        df['filename'] = name[start_index:end_index]
        
	# append dataframes to list
        df = df[final_cols]
        dataframes.append(df)

# concat everything in the list
end_df = pd.concat(dataframes).reset_index(drop=True)

# filter to only CDS
cds = end_df.loc[end_df['feature']=='CDS']

def sub_col(x, str_to_find, sep):
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

def split_filename(x, col_name):
    name_list = x.split('_')
    if col_name == 'mouse_name':
        return name_list[0]
    elif col_name == 'time_point':
        return name_list[1]
    elif col_name == 'isolate_number':
        return name_list[2]
    elif col_name == 'sample_name':
        return name_list[-4]
    elif col_name == 'sequencing_lane':
        return name_list[-3]
    elif col_name == 'bacteria':
        return name_list[-2]
    elif col_name == 'strain_type':
        return name_list[-1]

# insert attribute columns
attribute_cols = ['locus_tag', 'ID', 'inference', 'product',
                  'eC_number', 'name', 'db_xref', 'gene']
for col_name in attribute_cols:
        sub = col_name + '='
        cds.insert(loc=0, column=col_name,
                   value=end_df['attribute'].apply(sub_col,
                                                    str_to_find=sub,
                                                    sep=';'))
# insert filename columns
filename_cols = ['mouse_name', 'time_point', 'isolate_number',
                 'strain_type', 'bacteria', 'sequencing_lane',
                 'sample_name']
for colname in filename_cols:
        cds.insert(loc=0, column=colname,
                   value=cds['filename'].apply(split_filename,
                                               col_name=colname))

# get gene annotations
cds = cds.loc[cds['product']!='hypothetical protein']
gene_subset = cds.loc[cds['gene'].isin(cds['gene'].dropna().unique())]
cds_subset = cds.loc[cds['product'].str.contains('GN%3D')]
cds_subset['gene'] = cds_subset['product'].apply(sub_col, str_to_find='GN%3D', sep=' ')
annotated_genes = pd.concat([gene_subset, cds_subset])

# save df
annotated_genes.to_csv('annotated_genes.csv', index=False)
