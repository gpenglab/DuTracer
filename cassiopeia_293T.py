#This script provides basic workflow to extract the molecule_table and allele_table for further barcode analysis 
#and phylogenetic tree reconstruction
#Please check the cassiopeia pipeline for more information: https://cassiopeia-lineage.readthedocs.io/en/latest/user_guide.html
import pandas as pd
import numba
import cassiopeia as cas

import pandas as pd
import numpy as np
import numba
import pylab
import cassiopeia as cas
import matplotlib.pyplot as plt

%matplotlib inline

import seaborn as sns
sns.set_theme(style="ticks", rc = {"axes.spines.right": False, "axes.spines.top": False})

import sys
sys.path.append( '/data1/home/gdpeng/chengchen/dualproject/sample293T/' )
import utils

# The raw FASTQs
input_files = [
    "/data1/home/gdpeng/chengchen/data/293T/reads/293T_1.clean.fq.gz", "/data1/home/gdpeng/chengchen/data/293T/reads/293T_2.clean.fq.gz"
]

# The sample name, used for naming output files
name = 'sc293Tv2'
# Directory to output results, 
# the relative working path could not work in remote jupyter server running background, 
# only for local server running foreground
#output_directory = "../targetProcessData/293Tv2_preprocess_pipeline" 
#in remote jupyter server running background in VScode, use absolute path instead
output_directory = "/data1/home/gdpeng/chengchen/dualproject/targetProcessData/293Tv2_preprocess_pipeline"
# Path to the target site reference sequence in FASTA format
reference_filepath = "/data1/home/gdpeng/chengchen/dualproject/targetRawData/dual.ref.fasta"
#NOTE: THE PCT48.ref.fasta HAVE BEEN ADPATED TO OUR OWN DATA
# Number of threads to use, whenever parallelization is possible
n_threads = 8
# Whether to allow a single intBC to have multiple allele states
# For chemistries for which barcode == cell, this should be `False`.
allow_allele_conflicts = False
# Verbosity of logging
verbose = True

#cassiopeia.pp.setup(output_dir, verbose=verbose)
cas.pp.setup(output_directory, verbose=verbose)

bam_fp = cas.pp.convert_fastqs_to_unmapped_bam(
    input_files,
    chemistry='bgiv4',
    output_directory=output_directory,
    name=name,
    n_threads=n_threads
)

bam_fp = cas.pp.filter_bam(
    bam_fp,
    output_directory=output_directory,
    quality_threshold=10,
    n_threads=n_threads,
)

#home-made cellbc whitelist for cellbcs (cacontenated from barcode combinations from any two barcodes in the BGI barcode list), 
#cell barcode is the combination of two 10 sub-barcodes, please consult to BGI for the barcode list
#time-comsuming step
bam_fp = cas.pp.error_correct_cellbcs_to_whitelist(
    bam_fp,
    whitelist='/data1/home/gdpeng/chengchen/dualproject/targetScripts/BGIcellBCs.txt',
    output_directory=output_directory,
    n_threads=16,
)

umi_table = cas.pp.collapse_umis(
    bam_fp,
    output_directory=output_directory,
    max_hq_mismatches=3,
    max_indels=2,
    method='likelihood',
    n_threads=n_threads,
)

#need to convert the cell names first before filtering
#the following file was obtained after running BGI's DNBelab_C4 software
real_cells = pd.read_table("~/dualproject/targetRawData/sc293T_lineage_barcodeTranslate.txt",header=None)
real_cells.columns = ['cellBC','cell_ID']
umi_table = pd.merge(umi_table, real_cells, how='inner', on='cellBC')
umi_table['cellBC']=umi_table['cell_ID']

umi_table = cas.pp.resolve_umi_sequence(
    umi_table,
    output_directory=output_directory,
    min_umi_per_cell=10,
    min_avg_reads_per_umi=2.0,
    plot=False, #plot True face the issuse of parameter "basey"
)

umi_table = cas.pp.align_sequences(
    umi_table,
    ref_filepath=reference_filepath,
    gap_open_penalty=15,
    gap_extend_penalty=1,
    n_threads=n_threads,
)

umi_table = cas.pp.call_alleles(
    umi_table,
    ref_filepath=reference_filepath,
    barcode_interval=(24, 38),
    cutsite_locations=[96, 150, 182, 229],
    cutsite_width=12,
    context=True,
    context_size=5,
)

umi_table = umi_table[umi_table['intBC'].map(str).apply(len) == 14]

# tranfer cellID
cellBC2seq = real_cells.drop_duplicates('cell_ID').rename(columns={'cellBC':'cell_ID','cell_ID':"cellBC"})
umi_table = pd.merge(umi_table, cellBC2seq, how='inner', on='cellBC')

#collapse integration barcodes using starcode
STARCODE_TH = 2
#starcode collapse
(umi_table['cell_ID'] + umi_table['intBC']).to_csv("collapsing.txt",
                                                sep='\t',
                                                index=False,
                                                header=False)
import os
os.system("~/software/starcode/starcode -t 4 -d {} -s collapsing.txt > collapsing_result.txt".format(str(STARCODE_TH)))

final = pd.read_csv("collapsing_result.txt", sep='\t', header=None)
final['cellBC'] = final[0].apply(lambda x: x[:20])
final['intBC'] = final[0].apply(lambda x: x[20:])
final.rename(columns={1:"count"}, inplace = True)
final.drop(columns=[0], inplace = True)
os.system('rm collapsing_result.txt')
os.system('rm collapsing.txt')
umi_table = umi_table[(umi_table['cell_ID'] + "_" + umi_table['intBC']).isin(final['cellBC'] + "_" + final['intBC'])]

umi_table = cas.pp.error_correct_umis(
    umi_table,
    max_umi_distance=1,
    allow_allele_conflicts=allow_allele_conflicts,
    n_threads=n_threads,
)

umi_table = cas.pp.filter_molecule_table(
    umi_table,
    output_directory=output_directory,
    min_umi_per_cell=10, # same as above cutoff
    min_avg_reads_per_umi=2.0,
    min_reads_per_umi=max([np.percentile(umi_table['readCount'], 99) / 10, 5]), # minimum is 5. May tailed to your own datasets
    intbc_prop_thresh=0.5,
    intbc_umi_thresh=10,
    intbc_dist_thresh=1,
    doublet_threshold=0.35,
    allow_allele_conflicts=allow_allele_conflicts,
    plot=False,
)

allele_table = cas.pp.call_lineage_groups(
    umi_table,
    output_directory=output_directory,
    min_umi_per_cell=10,
    min_avg_reads_per_umi=2.0,
    min_cluster_prop=0.005,
    min_intbc_thresh=0.05,
    inter_doublet_threshold=0.35,
    kinship_thresh=0.25,
    plot=True,
)
