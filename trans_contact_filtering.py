#! /usr/bin/env python

import pandas as pd
import argparse


def run(args):
    in_filename = args.input
    out_filename = args.output
    threshold = args.threshold
    df = pd.read_csv(in_filename, delimiter='\t')
    df_trans = df.loc[df['chrom1'] != df['chrom2']]
    temp_list = pd.DataFrame()
    for index, row in df_trans.iterrows():
        if row[6] > threshold:
            temp_list = temp_list.append(row)
    new_list = temp_list[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count']]
    new_list.to_csv(out_filename, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="filter_for_high_trans_interaction_bins")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.add_argument("-thd", help="count threshold", dest="comp", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()