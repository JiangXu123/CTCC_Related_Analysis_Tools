#! /usr/bin/env python

import argparse
import tabix
import pandas as pd
import concurrent.futures



def run(args):
    input_pixel_file = args.input  #'/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/block0.tmp'
    bin_file = args.bin_file  # '/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/hg38_bined_chrom_200b.bed.gz'
    output_file = args.output
    tb = tabix.open(bin_file)  # load tabix with the modified bin file
    df = pd.read_csv(input_pixel_file, delimiter='\t')
    df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'contact_fq']
    df['bin1_id'] = None
    df['bin2_id'] = None
    for index, row in df.iterrows():    # tabix will return a iterator object, the content can only be obtained via iteration
        results1 = tb.query(row['chrom1'], row['start1'], row['end1'])
        results2 = tb.query(row['chrom2'], row['start2'], row['end2'])
        df.at[index, 'bin1_id'] = int(next(results1)[3])
        df.at[index, 'bin2_id'] = int(next(results2)[3])
    df.drop(columns=['chrom1','start1','end1','chrom2','start2','end2'])
    df.to_csv(output_file, sep='\t', header=False, index=False)


def main():
    parser = argparse.ArgumentParser(description="transform cooler dump -t --join pixel file to simple 'pin1_id','pin2_id',count form")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.add_argument("-mbin", help="modified bin file compressed with bgzip", dest="bin_file", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()