python convert_sth_to_fa.py 

fasta_dir=out_AZ20/msa_files/
end=PF*.fa
files_list=($fasta_dir$end)
new_end=.upper.fa
msa_dir=msa_files/
for f in "${files_list[@]}"
do
	name=$(basename $f .fa)
	awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' $f > $fasta_dir$name$new_end
done

python make_tree.py
