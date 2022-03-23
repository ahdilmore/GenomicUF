fasta_dir=fa_AZ51/
bed_dir=out_AZ51
genes_of_interest=($bed_dir/PF*)
end=_files/
end_fa='.fa'
for gene in "${genes_of_interest[@]}"
do
	output_path=$gene/
	file_array=($gene/*.bed)
	for file in "${file_array[@]}"
	do
		x=$(basename $file _AZ51_mut.bed)
		fasta_file=$fasta_dir$x*.fna
		bedtools getfasta -fo $output_path$x$end_fa -fi $fasta_file -bed $file -name
	done
	# merge all files for each gene to end .fa 
	echo $gene
	cat $output_path*.fa > $gene$end_fa
done
echo 'done!'
