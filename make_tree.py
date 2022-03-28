import qiime2 as q2
import glob
import os
from Bio import AlignIO 
from qiime2.plugins import phylogeny 


# for .sth file in the glob, do this 
sth_dir='out_AZ51/msa_files/'
str_files=glob.glob(sth_dir + '*.sth')
pfs_that_fail = ['PF04632.15.sth', 'PF05488.16.sth', 'PF05433.18.sth', 
		 'PF16085.8.sth', 'PF04830.16.sth', 'PF14375.9.sth', 
		 'PF12773.10.sth', 'PF17805.4.sth', 'PF11106.11.sth']

for f in str_files:
	basename = os.path.basename(f)
	print(basename)	
	if basename in pfs_that_fail: 
		continue
	outname = f.replace(".sth", ".fa")
	# file conversion from stockholm to fasta 
	AlignIO.convert(in_file=f, in_format="stockholm", 
			out_format="fasta", out_file=outname)

	# import aligned sequence into qiime 
	msa = q2.Artifact.import_data(type='FeatureData[AlignedSequence]',
				      view=outname)

	# make tree!!!
	basename = f.replace(".sth", "")
	tree = phylogeny.methods.fasttree(alignment=msa).tree
	q2.Artifact.export_data(tree, output_dir=basename+'_tree')
print('done!')
