import qiime2 as q2
import glob
from qiime2.plugins import phylogeny 


# for .sth file in the glob, do this 
fa_dir='out_AZ20/msa_files/'
fa_files=glob.glob(fa_dir + '*.upper.fa')

for f in fa_files:
	# import aligned sequence into qiime 
	msa = q2.Artifact.import_data(type='FeatureData[AlignedSequence]',
				      view=f)

	# make tree!!!
	basename = f.replace(".upper.fa", "")
	tree = phylogeny.methods.fasttree(alignment=msa).tree
	q2.Artifact.export_data(tree, output_dir=basename+'_tree')
print('done!')
