DATA_DIR=../../callaban/Projects/Adaptation/prokka_AZ20/
data_array=($DATA_DIR*L001)
OUTPUT_DIR=tRNA_files/
for FILE in "${data_array[@]}"
do
	file_name=${FILE#*$DATA_DIR}
	tRNA='_tRNA.gff'
	fa='_tRNA.fa'
	awk '$3=="tRNA"' $FILE/*.gff > $FILE/$file_name$tRNA
	bedtools getfasta -fi $FILE/*.fna -bed $FILE/$file_name$tRNA -fo $OUTPUT_DIR$file_name$fa
done
