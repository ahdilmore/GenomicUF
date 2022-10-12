import pandas as pd
import qiime2 as q2
import os
import glob
from pybedtools import BedTool
from Bio import AlignIO
from qiime2.plugins import phylogeny
BED_COLUMNS = ['seqname', 'start', 'end', 'name']

def _make_fasta_name(x, col_to_sort):
    return x[col_to_sort] + '_' + x['filename'] + '_' + str(x['rank'])

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

def _merge_features(wd, dir_to_merge):
    os.chdir(dir_to_merge)
    os.system("cat *.fa > merged.fa")
    os.chdir(wd)

def _extract_sequence(data_dict, feats_df, col_to_sort, out_path):
    """Function to output bed files of all annotated pfams/genes"""
    _prep_dataframe_for_bedtools(feats_df)

    wd = os.getcwd()
    for f_id in feats_df[col_to_sort].unique():
        f_df = feats_df.loc[feats_df[col_to_sort]==f_id]
        unique_files = feats_df['filename'].unique()

        if not os.path.exists(out_path + f_id):
            os.mkdir(out_path + f_id + '/')
        
        for f in unique_files:
            sub_df = f_df.loc[f_df['filename']==f]
            bed_file = BedTool.from_dataframe(sub_df[BED_COLUMNS])
            fasta_path = data_dict[f][1]
            bed_file.sequence(fi=fasta_path, name=True, fo=out_path+f_id+'/'+f+'.fa')
        
        _merge_features(wd, out_path+f_id)

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
            # set up variables for hmmfetch and hmmalign 
            dirname = os.path.dirname(path)
            pfid = os.path.basename(dirname)
            pf_hmm_out = dirname + '/hmmfile.hmm'
            msa_sth_out = dirname + '/msa.sth'
            # do hmmfetch, hmmalign 
            hmmfetch_cmd = "hmmfetch hmmer_assets/Pfam-A.hmm %s > %s"%(pfid, pf_hmm_out)
            hmmalign_cmd = "hmmalign -o %s --trim %s"%(msa_sth_out, path)
            os.system(hmmfetch_cmd)
            os.system(hmmalign_cmd)
            # convert sth to .fa file
            msa_fa_out = msa_sth_out.replace('.sth', '.mixedcase.fa')
            AlignIO.convert(in_file=msa_sth_out, in_format='stockholm', 
                            out_format='fasta', out_file=msa_fa_out)
            # convert all characters to uppercase for qiime2 import 
            upper_fa_out = msa_sth_out.replace('.sth', '.upper.fa')
            awk_input = "awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' %s > %s"%(msa_fa_out, upper_fa_out)
            os.system(awk_input)

def _tree_made(tree_path):
    return os.path.isfile(tree_path)

def _make_tree_aligned(path_to_aligned, out_directory):
    pfam_id = os.path.basename(os.path.dirname(path_to_aligned))
    if not _tree_made(out_directory + pfam_id + '/tree.nwk'): 
        # import aligned sequence into qiime2 & make tree 
        msa = q2.Artifact.import_data(type='FeatureData[AlignedSequence]', view=path_to_aligned)
        tree = phylogeny.methods.fasttree(alignment=msa).tree
        q2.Artifact.export_data(tree, out_directory + pfam_id)

def _make_tree_unaligned(path_to_unaligned, out_directory):
    feat_id = os.path.basename(os.path.dirname(path_to_unaligned))
    if not _tree_made(out_directory + feat_id + '/tree.nwk'):
        # import aligned sequence into qiime2 & make tree
        seq = q2.Artifact.import_data(type='FeatureData[Sequence]', view=path_to_unaligned)
        tree = phylogeny.pipelines.align_to_tree_mafft_fasttree(sequences=seq).tree
        q2.Artifact.export_data(tree, out_directory + feat_id)

def tree_construction(data_dict: dict, 
                      feats_df: pd.DataFrame,
                      out_path: str, 
                      construction_feature : str = 'Pfam') -> None:
    # process the sequences for tree construction 
    _process_sequences(data_dict, feats_df, out_path, construction_feature)
    # make subdirectory for tree information 
    if not os.path.exists(out_path + 'TreeData'):
            os.mkdir(out_path + 'TreeData')
    aligned_fastas = glob.glob(out_path + 'SequenceData/*/msa.upper.fa')
    unaligned_fastas = glob.glob(out_path + 'SequenceData/*/merged.fa')
    if aligned_fastas == []: 
        for path in unaligned_fastas: 
            _make_tree_unaligned(path, out_path+'TreeData/')
    else:
        for path in aligned_fastas: 
            _make_tree_aligned(path, out_path+'TreeData/')