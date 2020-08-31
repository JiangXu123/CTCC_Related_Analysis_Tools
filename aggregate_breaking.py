#! /usr/bin/env python

import pandas as pd
import time
import argparse


def run(args):
    start1 = time.perf_counter()
    input_file = args.input
    output_file = args.output
    df = pd.read_csv(input_file, delimiter='\t', names=['RNAME', 'start', 'breaking_pos', 'count'])
    df = df.groupby(['RNAME', 'start']).sum()['count'].reset_index()
    df['breaking_pos'] = df['start'] + 1
    df = df[['RNAME', 'start', 'breaking_pos', 'count']]
    df[['start', 'breaking_pos']] = df[['start', 'breaking_pos']].astype(int)
    df.to_csv(output_file, header=False, index=False, sep='\t')
    end1 = time.perf_counter()
    print(f'process finished in {round(end1 - start1, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="aggregate the counts")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()