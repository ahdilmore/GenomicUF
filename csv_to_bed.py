import pandas as pd

annotated_genes = pd.read_csv('annotated_genes.csv', 
			      dtype={'seqname':int, 'start':int,
				     'end':int})
genes_of_interest = ['garD']

for gene in genes_of_interest:
	df = annotated_genes.loc[annotated_genes['gene']==gene]
	unique_files = df['filename'].unique()

	for file in unique_files: 
		bed_name = gene + '_files/' + file + '.bed'
		sub_df = df.loc[df['filename']==file]
		sub_df = sub_df[['seqname', 'start', 'end']]
		sub_df.to_csv(bed_name, index=False,
			      sep='\t')
