import pandas as pd

annotated_genes = pd.read_csv('annotated_genes.csv', 
			      dtype={'seqname':int, 'start':int,
				     'end':int, 'db_xref':str,
				     'eC_number':str})
genes_of_interest = ['garD']
final_cols = ['seqname', 'source', 'name', 'start', 'end',
	      'score', 'strand', 'frame', 'attribute']

for gene in genes_of_interest:
	df = annotated_genes.loc[annotated_genes['gene']==gene]
	unique_files = df['filename'].unique()

	for file in unique_files: 
		bed_name = gene + '_files/' + file + '.bed'
		sub_df = df.loc[df['filename']==file]
		sub_df = sub_df[final_cols]
		sub_df.to_csv(bed_name, index=False,
			      header=False, sep='\t')
