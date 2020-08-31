#! /usr/bin/env python

import pandas as pd
import argparse
import time


def run(args):
    input_file = args.input
    output_file = args.output
    start_time = time.perf_counter()
    df = pd.read_csv(input_file, '\t')
    df.columns = ['RNAME', 'breaking_pos', 'end', 'count']
    df['start'] = df['end'] - 2
    df = df[['RNAME', 'start', 'breaking_pos', 'count']]
    df.to_csv(output_file, '\t', header=False, index=False)
    end_time = time.perf_counter()
    print(f'process finished in {round(end_time - start_time, 2)}')


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()