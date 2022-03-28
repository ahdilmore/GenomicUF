fasta_dir=fa_AZ51/
bed_dir=out_AZ51/
abbr=PF*
genes_of_interest=($bed_dir$abbr)
end_fa='.fa'
hmm_end=.hmm
sth_end=.sth
msa_dir=msa_files/
for gene in "${genes_of_interest[@]}"
do
	output_path=$gene/
	file_array=($gene/*.bed)
	for file in "${file_array[@]}"
	do
		x=$(basename $file _AZ51_mut.bed)
		fasta_file=$fasta_dir$x*.fna
		#bedtools getfasta -fo $output_path$x$end_fa -fi $fasta_file -bed $file -name
	done
	# merge all files for each gene to end .fa 
	echo $gene
	cat $output_path*.fa > $gene$end_fa
	# generate hmm file for that pfam and perform msa
	pfam=$(basename $gene .fa)
	hmmfetch Pfam-A.hmm $pfam > $bed_dir$msa_dir$pfam$hmm_end 
	hmmalign -o $bed_dir$msa_dir$pfam$sth_end --trim $bed_dir$msa_dir$pfam$hmm_end $gene$end_fa
done
echo 'done!'
