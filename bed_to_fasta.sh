fasta_dir=../../callaban/Projects/Adaptation/prokka_AZ20/
genes_of_interest=(garD)
end=_files/
end_fa='.fa'
for gene in "${genes_of_interest[@]}"
do
	output_path=$gene$end
	file_array=($output_path*.bed)
	
	for file in "${file_array[@]}"
	do
		name="`awk -F'_files/|_Nissle' '{print $2}' <<< "$file"`"
		fasta_file=$fasta_dir$name/*.fna
		echo $fasta_file
		echo $file
		bedtools getfasta -fo $output_path$name$end_fa -fi $fasta_file -bed $file
	done
	# merge all files for each gene to end .fa 
	cat $output_path*.fa > $gene$end_fa
done
