import glob
import os
from Bio import AlignIO 

sth_dir='out_AZ20/msa_files/'
str_files=glob.glob(sth_dir + '*.sth')

for f in str_files:
	basename = os.path.basename(f)
	outname = f.replace(".sth", ".fa")
	# file conversion from stockholm to fasta 
	AlignIO.convert(in_file=f, in_format="stockholm", 
			out_format="fasta", out_file=outname)
