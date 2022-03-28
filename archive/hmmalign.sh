out_dir=out_AZ51/
fasta_dir=($out_dir*.fa)
hmm_end=.hmm
sth_end=.sth
msa_dir=msa_files/
for pfam in "${fasta_dir[@]}"
do
	fn=$(basename $pfam .fa)
	echo $fn
	hmmfetch Pfam-A.hmm $fn > $out_dir$msa_dir$fn$hmm_end
	hmmalign -o $out_dir$msa_dir$fn$sth_end --trim $out_dir$msa_dir$fn$hmm_end $pfam 
done
