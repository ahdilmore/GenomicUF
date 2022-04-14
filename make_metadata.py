import pandas as pd
import numpy as np 

def split_filename(x, col_name):
	name_list = x.split('_')
	if col_name == 'mouse_name':
		return name_list[0]
	elif col_name == 'time_point':
		return name_list[1]
	elif col_name == 'isolate_number':
		return name_list[2]
	elif col_name == 'additional_name':
		if(len(name_list)==7):
			return name_list[-4]
		return np.nan
	elif col_name == 'sequencing_lane':
		return name_list[-3]
	elif col_name == 'bacteria':
		return name_list[-2]
	elif col_name == 'strain_type':
		return name_list[-1]

def make_dataframe(sample_names, end_cols):
	df = pd.DataFrame(data={'sample_name': sample_names})
	for colname in end_cols:
		df.insert(loc=0, column=colname,
			  value=df['sample_name'].apply(split_filename,
							col_name=colname))
	df = df.set_index('sample_name')
	return df

filename_cols = ['mouse_name', 'time_point', 'isolate_number',
                 'strain_type', 'bacteria', 'sequencing_lane',
                 'additional_name']
start_df = pd.read_csv('pfam_annots_AZ51.csv')
end_df = make_dataframe(start_df['filename'].unique(), filename_cols)
end_df.to_csv('AZ51_metadata.csv')
