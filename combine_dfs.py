import pandas as pd 
import glob 

files = glob.glob('gff_files/*.gff')
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
	start_index = f.index('/') + 1
	end_index = f.index('.')
	df['filename'] = f[start_index:end_index]
	
	# append dataframes to list
	df = df[final_cols]
	dataframes.append(df)	

# concat everything in the list
df = pd.concat(dataframes).reset_index(drop=True)

# remove blank lines between each concat
df = df.dropna()

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
        return name_list[3]
    elif col_name == 'sequencing_lane':
        return name_list[4]
    elif col_name == 'bacteria':
        return name_list[5]
    elif col_name == 'strain_type':
        return name_list[6]

# insert attribute columns
df.insert(loc=0, column='locus_tag', value=df['attribute'].apply(sub_col, str_to_find='locus_tag=', sep=';'))
df.insert(loc=0, column='ID', value=df['attribute'].apply(sub_col, str_to_find='ID=', sep=';'))
df.insert(loc=0, column='inference', value=df['attribute'].apply(sub_col, str_to_find='inference=', sep=';'))
df.insert(loc=0, column='product', value=df['attribute'].apply(sub_col, str_to_find='product=', sep=';'))
df.insert(loc=0, column='eC_number', value=df['attribute'].apply(sub_col, str_to_find='eC_number=', sep=';'))
df.insert(loc=0, column='name', value=df['attribute'].apply(sub_col, str_to_find='Name=', sep=';'))
df.insert(loc=0, column='db_xref', value=df['attribute'].apply(sub_col, str_to_find='db_xref=', sep=';'))
df.insert(loc=0, column='gene', value=df['attribute'].apply(sub_col, str_to_find='gene=', sep=';'))

# insert filename columns
df.insert(loc=0, column='mouse_name', value=df['filename'].apply(split_filename, col_name='mouse_name'))
df.insert(loc=0, column='time_point', value=df['filename'].apply(split_filename, col_name='time_point'))
df.insert(loc=0, column='isolate_number', value=df['filename'].apply(split_filename, col_name='isolate_number'))
df.insert(loc=0, column='strain_type', value=df['filename'].apply(split_filename, col_name='strain_type'))
df.insert(loc=0, column='bacteria', value=df['filename'].apply(split_filename, col_name='bacteria'))
df.insert(loc=0, column='sequencing_lane', value=df['filename'].apply(split_filename, col_name='sequencing_lane'))
df.insert(loc=0, column='sample_name', value=df['filename'].apply(split_filename, col_name='sample_name'))

# filter to only CDS & remove hypothetical proteins
cds = df.loc[(df['feature']=='CDS') & (df['product']!='hypothetical protein')]

# save df
cds.to_csv('cds_intermediate.csv')
