#! /usr/bin/env python
# this program convert the float values in previous program to int
import pandas as pd
import time
import argparse


def run(args):
    start1 = time.perf_counter()
    input_file = args.input_f
    output_file = args.output_f
    df = pd.read_csv(input_file, delimiter='\t', names=['RNAME', 'start', 'end', 'counts'])
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    df.to_csv(output_file, header=False, index=False, sep='\t')
    end1 = time.perf_counter()
    print(f'process finished in {round(end1 - start1, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="generate the bed file for breaking analysis for chromosome 1")
    parser.add_argument("-in", help="input file", dest="input_f", type=str, required=True)
    parser.add_argument("-out", help="output file", dest="output_f", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()