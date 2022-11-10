import pandas as pd
import numpy as np
#import qiime2 as q2
import os
import subprocess
import glob
import skbio
from pybedtools import BedTool
from Bio import AlignIO
#from qiime2.plugins import phylogeny
BED_COLUMNS = ['seqname', 'start', 'end', 'name']

def run_command(cmd, output_fp=None):
    if output_fp is None: 
        subprocess.run(cmd)
    else: 
        with open(output_fp, 'w') as f: 
            subprocess.run(cmd, stdout=f)

def _make_fasta_name(x, col_to_sort):
    return x[col_to_sort] + '_' + x['filename'] + '_' + str(x['rank'])

def _file_made(path):
    return os.path.isfile(path)

def _prep_dataframe_for_bedtools(feats_df, col_to_sort): 
    # ensure that feats_df.start and feats_df.end are integers
    feats_df.start = feats_df.start.astype('int')
    feats_df.end = feats_df.end.astype('int')
    # get rank for placing multiple sequences in tree
    feats_grouped = feats_df.groupby([col_to_sort, 'filename'])
    feats_df.insert(loc=0, column='rank', 
                    value=feats_grouped['start'].rank().astype('int'))
    # make name of fasta file using ranks calculated above 
    fasta_names = feats_df.apply(_make_fasta_name, args=(col_to_sort, ), axis=1)
    feats_df.insert(loc=0, column='name', value=fasta_names)

def _merge_features(dir_to_merge):
    files_to_merge = glob.glob(dir_to_merge + '*.fa')
    run_command(["cat"] + files_to_merge, dir_to_merge + 'merged.fa')

def _extract_sequence(data_dict, feats_df, col_to_sort, out_path):
    """Function to output bed files of all annotated pfams/genes"""
    _prep_dataframe_for_bedtools(feats_df, col_to_sort)

    for f_id in feats_df[col_to_sort].unique():
        f_df = feats_df.loc[feats_df[col_to_sort]==f_id]
        unique_files = feats_df['filename'].unique()

        if not os.path.exists(out_path + f_id):
            os.mkdir(out_path + f_id)
        
        if not _file_made(out_path+f_id+'/merged.fa'):
            for f in unique_files:
                if not _file_made(out_path+f_id+'/'+f+'.fa'):
                    sub_df = f_df.loc[f_df['filename']==f]
                    bed_file = BedTool.from_dataframe(sub_df[BED_COLUMNS])
                    fasta_path = data_dict[f][1]
                    bed_file.sequence(fi=fasta_path, name=True, fo=out_path+f_id+'/'+f+'.fa')
            _merge_features(out_path+f_id+'/')
        

def _hmmer_alignment(path_to_merged):
    # set up variables for hmmfetch and hmmalign 
    dirname = os.path.dirname(path_to_merged)
    pfid = os.path.basename(dirname)
    pf_hmm_out = dirname + '/hmmfile.hmm'
    msa_sth_out = dirname + '/msa.sth'
    # if files don't exist already, do hmmfetch, hmmalign 
    if not _file_made(pf_hmm_out):
        hmmfetch_cmd = ["hmmfetch", "../../preprocessing/hmmer_assets/Pfam-A.hmm", pfid]
        run_command(hmmfetch_cmd, pf_hmm_out)
    if not _file_made(msa_sth_out):
        hmmalign_cmd = ["hmmalign", "-o", msa_sth_out, "--trim", pf_hmm_out, path_to_merged]
        run_command(hmmalign_cmd)
        # convert sth to .fa file
        msa_fa_out = msa_sth_out.replace('.sth', '.mixedcase.fa')
        AlignIO.convert(in_file=msa_sth_out, in_format='stockholm', 
                        out_format='fasta', out_file=msa_fa_out)
        # convert all characters to uppercase for qiime2 import 
        awk_cmd = ["awk", 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}', msa_fa_out]
        run_command(awk_cmd, msa_sth_out.replace('.sth', '.upper.fa'))

def _process_sequences(data_dict, feats_df, out_path, col_to_sort):
    # make subdirectory: sequence data 
    if not os.path.exists(out_path + 'SequenceData'):
            os.mkdir(out_path + 'SequenceData')
    # extract all the subsequences, merge to file for msa
    _extract_sequence(data_dict, feats_df, col_to_sort, out_path + 'SequenceData/')
    
    if col_to_sort == 'Pfam':
        # do HMM Alignment  
        merged_paths = glob.glob(out_path + 'SequenceData/*/merged.fa') 
        for path in merged_paths:
            _hmmer_alignment(path)

def _mafft(alignment_fp, n_threads=1):
    cmd = ['mafft', '--preservecase', '--inputorder', '--thread', str(n_threads)]
    run_command(cmd, alignment_fp)

def _fasttree(aligned_fp, tree_fp):
    cmd = ['FastTree', '-quote', '-nt', aligned_fp]
    run_command(cmd, tree_fp)


def _make_biom_table(tree_path, gene_name):
    tree = skbio.io.read(tree_path, format="newick", into=skbio.TreeNode)
    # intialize lists of nodes and sample names 
    node_names = []
    sample_names = []
    # make list of node names 
    for node in tree.tips():
        node_names.append(node.name)
    # make list of the sample names 
    for name in node_names:
        final_index = name.index("::")
        sample_name = name[len(gene_name)+1:final_index]
        if sample_name not in sample_names:
            sample_names.append(sample_name)
    # put array together 
    empty_array = np.zeros((len(sample_names), len(node_names)))
    for i in range(len(sample_names)):
        for j in range(len(node_names)):
            if sample_names[i] in node_names[j]:
                empty_array[i][j] = 1
    biom_table = pd.DataFrame(data=empty_array, columns=node_names,
                              index=sample_names)
    biom_table = biom_table.T
    # save biom table
    txt_path = tree_path.replace('tree.nwk', 'table.txt')
    biom_table.to_csv(txt_path, sep='\t')
    # convert to .biom format for later use 
    biom_path = tree_path.replace('tree.nwk', 'table.biom')
    convert_txt_cmd = ["biom", "convert", "-i", txt_path, "-o", biom_path, "--to-hdf5"]
    run_command(convert_txt_cmd)

def _make_tree_aligned(path_to_aligned, out_directory):
    pfam_id = os.path.basename(os.path.dirname(path_to_aligned))
    # check that the directory has been made
    if not os.path.exists(out_directory + pfam_id):
            os.mkdir(out_directory + pfam_id)
    # if tree has not been made, run fasttree with the MSA 
    if not _file_made(out_directory + pfam_id + '/tree.nwk'): 
        _fasttree(path_to_aligned, out_directory + pfam_id + '/tree.nwk')
    # make biom table based on the constructed tree
    if not _file_made(out_directory + pfam_id + '/table.biom'):
        _make_biom_table(out_directory + pfam_id + '/tree.nwk', pfam_id)

'''
def _make_tree_unaligned(path_to_unaligned, out_directory):
    feat_id = os.path.basename(os.path.dirname(path_to_unaligned))
    if not _file_made(out_directory + feat_id + '/tree.nwk'):
        # perform sequence alignment with mafft
        seq = q2.Artifact.import_data(type='FeatureData[Sequence]', view=path_to_unaligned)
        tree = phylogeny.pipelines.align_to_tree_mafft_fasttree(sequences=seq).tree
        q2.Artifact.export_data(tree, out_directory + feat_id)
    # make biom table based on the constructed tree
    if not _file_made(out_directory + feat_id + '/table.biom'):
        _make_biom_table(out_directory + feat_id + '/tree.nwk', feat_id)
 '''

def tree_construction(data_dict: dict, 
                      feats_df: pd.DataFrame,
                      out_path: str,
                      construction_feature : str = 'Pfam') -> None:
    '''This pipeline constructs per-gene trees that will be used downstream
    in single UniFrac or meta UniFrac calculations. 
    data_dict: dictionary where keys are sample names and values are a tuple of the
    .gff filepath and the .fa filepath for that sample
    feats_df: processed pd.DataFrame with filtered annotations that was made in preprocessing steps
    out_path: an output directory where the subsetted sequencing data and the tree data will be written.
    construction_feature: determines whether we will be using Pfam or gene to construct trees.'''
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # check that bed columns & construction_feature are in the feats_df
    check_cols = [col for col in BED_COLUMNS if col!='name'] + [construction_feature]
    for col in check_cols:
        if col not in feats_df.columns: 
            raise ValueError('The ' + col + ' column is not present in input dataframe.')
    # process the sequences for tree construction 
    _process_sequences(data_dict, feats_df, out_path, construction_feature)
 
    # make subdirectory for tree information 
    if not os.path.exists(out_path + 'TreeData'):
            os.mkdir(out_path + 'TreeData')
    aligned_fastas = glob.glob(out_path + 'SequenceData/*/msa.upper.fa')
    #unaligned_fastas = glob.glob(out_path + 'SequenceData/*/merged.fa')
    #if aligned_fastas == []: 
    #    for path in unaligned_fastas: 
    #       _make_tree_unaligned(path, out_path+'TreeData/')
    #else:
    for path in aligned_fastas: 
        _make_tree_aligned(path, out_path+'TreeData/')