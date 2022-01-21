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
end_df = pd.concat(dataframes).reset_index(drop=True)
end_df.to_csv('very_big_df.csv')
