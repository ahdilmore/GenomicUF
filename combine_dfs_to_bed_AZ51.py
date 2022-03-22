import pandas as pd
import numpy as np
import glob
import os

# SET PARAMS: 
files_dir = '../../../panfs/callaban/prokka_AZ51/'
files_pattern = '*S*/*.gff'
files_cols = ['filename', 'seqname', 'source', 'feature', 'start',
              'end', 'score', 'strand', 'frame', 'attribute']
attribute_cols = ['ID', 'Parent', 'eC_number', 'Name',
                  'dbxref', 'gene', 'inference', 'locus_tag',
                  'product', 'protein_id']
final_cols = ['seqname', 'start', 'end', 'name']
to_remove = ['../../../panfs/callaban/prokka_AZ51/46R_Month1_13_S1816/46R_Month1_13_S1816_AZ51_mut.gff']
out_path = 'out_AZ51/'

def concat_annotations(files_dir, files_pattern, files_cols, to_remove): 
	"""This function takes directory containing several 
	.gff/.csv files and concatenates them into a pd.DataFrame
	to make parsing easier for downstream steps."""
	files = glob.glob(files_dir + files_pattern)
	
	# remove files that are causing issues
	for i in to_remove: 
		files.remove(i)
	
	# convert each gff to dataframe; make list of dfs 
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
				  data = array).dropna()
		name = f[len(files_dir):]
		start_index = name.index('/') + 1
		end_index = name.index('.')
		df.insert(loc=0, column='filename',
			  value=name[start_index:end_index])
        
        	# append dataframes to list
		df = df[files_cols]
		dataframes.append(df)
    
    	# concat everything in the list
	return pd.concat(dataframes).reset_index(drop=True)

def sub_col(x, str_to_find, sep):
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
    
def split_attribute(df, cols_to_insert):
	"""Helper function to split the attribute column in gff 
	files so that more meaningful information can be extracted"""
	# insert columns 
	for col_name in cols_to_insert:
		sub = col_name + '='
		df.insert(loc=0, column=col_name,
			  value=df['attribute'].apply(sub_col,
						      str_to_find=sub,
						      sep=';'))
	# remove attribute
	df.drop(columns=['attribute'], inplace=True)
	return df

def extract_pfam(files_dir, files_pattern, files_cols, to_remove, 
		 attribute_cols):
	"""Function to perform all steps of analysis for you. Takes in a
	list of gff / csv files and outputs a pandas dataframe containing 
	the Pfam annotations."""
    	# concat dataframes
	concat = concat_annotations(files_dir, files_pattern, files_cols,
				    to_remove)
    
   	# filter to CDS 
	cds = concat.loc[concat['feature']=='CDS']
    
    	# split the attribute column
	cds = split_attribute(cds, attribute_cols)
    
    	# pfam filter
	pfam = cds.loc[cds['inference'].str.contains('Pfam')]
	pfam.insert(loc=0, column='Pfam',
		    value=pfam['inference'].apply(sub_col,
						  str_to_find='Pfam-A:', 
						  sep=':'))
	return pfam

def get_id(x, str_to_find):
	if str_to_find in x:
		end_index = x.index(str_to_find)
		return x[:end_index]
	return x

def make_fasta_name(x):
	return x['pfam_id'] + '_' + x['filename'] + '_' + str(x['rank'])

def make_bed_files(pfam, final_cols, out_path):
	"""Function to output bed files of all pfam annotated-genes"""
	pfam.insert(loc=0, column='rank', value=pfam.groupby(['Pfam', 'filename'])['start'].rank().astype('int'))
	pfam.insert(loc=0, column='name', value=pfam.apply(make_fasta_name, axis=1))

	for pfam_id in pfam['Pfam'].unique():
		df = pfam.loc[pfam['Pfam']==pfam_id]
		unique_files = df['filename'].unique()
		if not os.path.exists(out_path + pfam_id):
			os.mkdir(out_path + pfam_id + '/')
		for f in unique_files:
			bed_name = out_path + pfam_id + '/' + f + '.bed'
			sub_df = df.loc[df['filename']==f]
			sub_df[final_cols].to_csv(bed_name, index=False, header=False, sep='\t')


# combine dfs
pfam = extract_pfam(files_dir, files_pattern, files_cols, to_remove, attribute_cols)
pfam.to_csv('pfam_annots_AZ51.csv', index=False)

# make the bed files
pfam = pd.read_csv('pfam_annots_AZ51.csv')
make_bed_files(pfam, final_cols, out_path) 
