python make_biom.py

table_dir=out_AZ20/msa_files/*tree/
table_ext=*.txt
out_ext=/table.biom
tables=($table_dir$table_ext)
for table in "${tables[@]}"
do
	out_name=$(dirname $table)
	echo $out_name$out_ext
	biom convert -i $table -o $out_name$out_ext --to-hdf5
done
