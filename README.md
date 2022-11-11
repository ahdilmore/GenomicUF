# GenomicUF
Utilizes microbial genomic annotations - such as those output from [prokka](https://github.com/tseemann/prokka) - to construct per-gene trees, which may then be used to perform Unweighted UniFrac calculations on a per-gene basis. We can perform either a single-gene analysis, yielding an effect size for each unique gene or a multi-gene analysis i.e. [Meta UniFrac](https://github.com/biocore/unifrac/blob/master/unifrac/_meta.py) yielding an effect size for several groups of genes. Note that the latter does not scale well computationally and only small numbers of combinations are recommended at this point.  

## Installation 
```
conda install -c bioconda bedtools hmmer fasttree unifrac -c anaconda scikit-bio
git clone https://github.com/ahdilmore/GenomicUF.git
cd GenomicUF 
pip install -e .
```
