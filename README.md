# GenomicUF
Utilizes microbial genomic annotations - such as those output from [prokka](https://github.com/tseemann/prokka) - to construct per-gene trees, which may then be used to perform Unweighted UniFrac calculations on a per-gene basis. We can perform either a single-gene analysis, yielding an effect size for each unique gene or a multi-gene analysis i.e. [Meta UniFrac](https://github.com/biocore/unifrac/blob/master/unifrac/_meta.py) yielding an effect size for several groups of genes. Note that the latter does not scale well computationally and only small numbers of combinations are recommended at this point.  

## Installation 
```
conda install -c bioconda bedtools hmmer fasttree unifrac -c anaconda scikit-bio
git clone https://github.com/ahdilmore/GenomicUF.git
cd GenomicUF 
pip install -e .
```

## What is Genomic UniFrac?
Genomic UniFrac implements the [UniFrac algorithm](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1317376/) on whole genome sequences. More than one microbial species may be sequenced (i.e. the data may be metagenomic in neature), but whole genome sequences are available. This is distinct from Meta UniFrac (or Metagenomic UniFrac) which operates on sequencing data obtained from a mixed microbial community but whole genome sequences are not necessarily inferred. 

## Workflow: 

The input to Genomic UniFrac is either a data dictionary or a glob pattern. 

* The keys to the data dictionary are the sample names. Each value is a tuple. The first element in the tuple is the path to the .gff file containing annotations and the second element in the tuple is the path to the .fa file containing genomic sequence information. 
* Alternatively, one can specify a glob pattern to fetch the .gff and .fa files in a given directory. In our example case, the name of each directory is the name of the sample for which it contains annotation and genomic information. 

**Preprocessing steps** 
* Inputs to the Genomic UniFrac pipeline go through the **process_data_dict** function, which: 
  * Checks that only one of the data dictionary or glob pattern is present, not both
  * If a glob pattern is given, checks that one .gff file and one .fa file are present in each directory - and raises an exception if that is not the case. This step also converts the glob pattern into a data dictionary for downstream processing.  
  * If a data dictionary is given, checks that the data dictionary is in the correct format and checks that the filepaths given in the data dictionary exists. 
* The annotations files in the data dictionary are concatenated in the **concat_annotations** function. The idea behind this function is that the output dataframe can be sorted for elements of interest, which will be run through the Genomic UniFrac analysis pipeline. In this step, we also check that the .gff files actually contain genomic annotation information (a warning is raised and the sample is skipped in this case) and we check that the .gff files provided are of the correct format (raises an exception if not)

  Some example filtering perfomed after concatenating the annotations together include: 
  * The **filter_annotations** function, which allows one to filter to a specific feature of interest. The default behavior is to filter to coding sequences only (CDS), but other relevant examples include tRNA, etc. 
  * The **extract_pfam** function identifies entries that are flagged with their Pfam domain.
* The **wrapper_function** does all these preprocessing steps together. 

**Tree construction steps** 
Inputs to tree construction are the data dictionary and feature dataframe created in the preprocessing steps, as well as a metadata file for the samples of interest. We provide the function with the feature we want to construct the trees based on (e.g. Pfam trees or gene trees). 
 * The **_extract_sequence** function writes out the sequence infomation for each entry in an input GFF annotations table using the BedTools module. 
 * The **_process_sequence** function 1) makes an output directory for sequence information for each gene in the annotation file, 2) iterates through each gene and each sample in the annotation file and runs _extract_sequence for each, and 3) concatenates all the sequence information together into one .fa file for a multiple sequence alignment. If Pfams are the feature of interest here, the HMMER multiple sequence alignment is done within this step. 
 * If Pfams are the feature of interest, we run **_make_tree_aligned**, which imports the HMMER multiple sequence alignment into QIIME2 Phylogeny and runs MAFFT fasttree. If genes are the feature of interest, we run **_make_tree_unaligned** which imports the merged sequences from **_process_sequence** into QIIME2 Phylogeny and runs the Align to Tree MAFFT fasttree pipeline. Within each of these tree construction steps, after the tree is constructed, we call **_make_biom_table** which converts the information in the tree into a presence/absence biom table so that we can run UniFrac. 

**Analysis steps** 
* The **single_genes** function runs UniFrac and generates an effect size for each individual gene/Pfam that has previously been preprocessed. The col_of_interest parameter is the variable we will be generating an effect size for. 
