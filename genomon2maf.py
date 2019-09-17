#!/usr/bin/env python
"""
Generate an input MAF file for MutSigCV by converting SNV/INDEL calls by the Genomon2 pipeline.
The input is a multi-sample variant call in post_analysis/ directory.
Usage: python genomon2maf.py merge_mutation_filt.txt
"""
__author__ = "Masashi Fujita <mssfjt@gmail.com>"
__version__ = "0.0.1"
__date__ = "September 17, 2019"

import sys
import argparse
import pandas as pd
import numpy as np

func_dict = {
    'intronic': 'Intron',
    'ncRNA_exonic': 'RNA',
    'ncRNA_intronic': 'RNA',
    'ncRNA_splicing': 'RNA',
    'splicing': 'Splice_Site',
    'UTR3': '3\'UTR',
    'UTR5': '5￿\'UTR'
}
exonic_func_dict = {
    'frameshift deletion': 'Frame_Shift_Del',
    'frameshift insertion': 'Frame_Shift_Ins',
    'nonframeshift deletion': 'In_Frame_Del',
    'nonframeshift insertion': 'In_Frame_Ins',
    'nonsynonymous SNV': 'Missense_Mutation',
    'stopgain': 'Nonsense_Mutation',
    'stoploss': 'Nonstop_Mutation',
    'synonymous SNV': 'Silent'
}


def fetch_variant_classification(row):
    func = row['Func.refGene']
    exonicFunc = row['ExonicFunc.refGene']
    if func in 'exonic' and exonicFunc in exonic_func_dict:
        return exonic_func_dict[exonicFunc]
    elif func in func_dict:
        return func_dict[func]
    else:
        return np.nan


#
# parse args
#
parser = argparse.ArgumentParser(description='Convert Genomon2 SNV/INDEL calls to MAF.')
parser.add_argument('infile', metavar='genomon_file', help='Genomon2 SNV/INDEL file for a single sample')
parser.add_argument('--out', '-o', metavar='MAF', help='Output MAF file [default: stdout]')
args = parser.parse_args()

#
# open input
#
df = pd.read_csv(args.infile, sep='\t', skiprows=3, dtype={'Chr': 'str'}, low_memory=False)

#
# format into MAF
#
df.rename(columns={
    'Gene.refGene': 'Hugo_Symbol',
    'Chr': 'Chromosome',
    'Start': 'Start_position',
    'End': 'End_position',
    'Ref': 'Reference_Allele',
    'Alt': 'Tumor_Seq_Allele1',
    'id': 'Tumor_Sample_Barcode'
    }, inplace=True)
df['NCBI_Build'] = '37'
df['Tumor_Seq_Allele2'] = df['Tumor_Seq_Allele1']
df['Variant_Classification'] = df.apply(lambda row: fetch_variant_classification(row), axis=1)
df = df[['Hugo_Symbol', 'NCBI_Build', 'Chromosome', 'Start_position', 'End_position', 'Variant_Classification',
         'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'Tumor_Sample_Barcode',
         'ExonicFunc.refGene', 'Func.refGene']]

#
# save
#
if args.out is not None:
    f_out = open(args.out, 'w')
else:
    f_out = sys.stdout
df.to_csv(f_out, sep='\t', header=True, index=False, na_rep="NA")
f_out.close()
