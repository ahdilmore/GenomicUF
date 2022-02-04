# GenomicUF
Scripts to implement Genomic UniFrac for Celeste's adaptation project. The goal is to eventually generalize this codebase so that Genomic UniFrac can more readily be applied to many studies. 

## What is Genomic UniFrac?
Genomic UniFrac implements the [UniFrac algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1317376/) on whole genome sequences. More than one microbial species may be sequenced (i.e. the data may be metagenomic in neature), but whole genome sequences are available. This is distinct from Meta UniFrac (or Metagenomic UniFrac) which operates on sequencing data obtained from a mixed microbial community but whole genome sequences are not necessarily inferred. 

## Required inputs
* .fasta / .fa / .fna files containing genome sequences 
* .gff (general feature format; tab delimited file containing genomic features) files for each genome sequence

## General analysis pipeline
* Perform multiple sequence alignment for sequences with same annotations across different samples/timepoints
  * Quality control on GFF files obtained from sequences
  * Concatenate GFF files into a pandas DataFrame for ease of parsing information
  * Slice dataframe by gene of interest
  * Extract individual genes from genome fasta files; concat those from every genome file into one fasta file 
  * Align genomic sequences using [HMMER3](http://hmmer.org/)
* Construct per-gene phylogenetic trees with [MAFFT fasttree](https://docs.qiime2.org/2021.11/plugins/available/phylogeny/align-to-tree-mafft-fasttree/?highlight=mafft%20fast%20tree)  
* Calculate Meta UniFrac [code](https://github.com/biocore/unifrac/blob/077fca46bd)
  * Hypothetically, there will be some numerical optimization built in to handle the large number of combinations of genes 

## Current required steps
1) Run combine_dfs.py to generate general dataframe for searching genes 
2) Run csv_to_bed.py to generate bed files for particular genes of interest 
3) Run bed_to_fasta.sh to slice sequencing files for input to HMMER alignment 

## HAZEL TODO
1) Troubleshoot the bed files that do not have valid .fna file associated 
 * I think this is mostly solved - will probably have to play around with the order of the columns, but besides that I think everything else should work
2) Try HMMER alignment with one gene of interest 
3) Build gene tree for one gene of interest
