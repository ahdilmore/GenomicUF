import pandas as pd
import numpy as np
import glob
import os

# SET PARAMS: 
files_dir = '../../../panfs/callaban/prokka_AZ20/'
files_pattern = '*S*/*.gff'
files_cols = ['filename', 'seqname', 'source', 'feature', 'start',
              'end', 'score', 'strand', 'frame', 'attribute']
attribute_cols = ['ID', 'Parent', 'eC_number', 'Name',
                  'dbxref', 'gene', 'inference', 'locus_tag',
                  'product', 'protein_id']
final_cols = ['seqname', 'start', 'end', 'name']
to_remove = []
out_path = 'out_AZ20/'

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

def make_fasta_name(x, col_to_sort):
	return x[col_to_sort] + '_' + x['filename'] + '_' + str(x['rank'])

def make_bed_files(pfam, final_cols, out_path, col_to_sort):
	"""Function to output bed files of all annotated pfams/genes"""
	pfam.insert(loc=0, column='rank', value=pfam.groupby([col_to_sort, 'filename'])['start'].rank().astype('int'))
	pfam.insert(loc=0, column='name', value=pfam.apply(make_fasta_name, args=(col_to_sort, ), axis=1))

	for pfam_id in pfam[col_to_sort].unique():
		df = pfam.loc[pfam[col_to_sort]==pfam_id]
		unique_files = df['filename'].unique()
		print(pfam_id)
		if not os.path.exists(out_path + pfam_id):
			os.mkdir(out_path + pfam_id + '/')
		for f in unique_files:
			bed_name = out_path + pfam_id + '/' + f + '.bed'
			sub_df = df.loc[df['filename']==f]
			sub_df[final_cols].to_csv(bed_name, index=False, header=False, sep='\t')


# combine dfs
#pfam = extract_pfam(files_dir, files_pattern, files_cols, to_remove, attribute_cols)
#pfam.to_csv('pfam_annots_AZ20.csv', index=False)

# make the bed files
#pfam = pd.read_csv('pfam_annots_AZ20.csv')
#pfam_over_5 = pfam["Pfam"].value_counts()[pfam["Pfam"].value_counts() > 5].index
#pfam = pfam.loc[pfam["Pfam"].isin(pfam_over_5)]
#make_bed_files(pfam, final_cols, out_path)

# FOR GENE COMPARISON: 
#full_df = concat_annotations(files_dir, files_pattern, files_cols, to_remove)
#genes_df = full_df.loc[full_df["feature"] == "CDS"]
#genes_df = split_attribute(genes_df, attribute_cols)
#genes_df = genes_df.loc[genes_df["gene"].notna()]
#genes_df.to_csv('AZ20_full_df.csv', index=False)

genes_df = pd.read_csv("AZ20_full_df.csv")
# genes present must be present in over 100 isolates (~10%)
genes_over_100 = genes_df["gene"].value_counts()[genes_df["gene"].value_counts() > 100].index
genes_df = genes_df.loc[genes_df["gene"].isin(genes_over_100)] 
out_path_new = 'genes_out_AZ20_filt/'
make_bed_files(genes_df, final_cols, out_path_new, "gene") 
 
