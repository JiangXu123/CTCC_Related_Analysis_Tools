#! /usr/bin/env python

import pandas as pd
import time
import concurrent.futures
import argparse
# this code was writen in the hope of finding the breaking_point using multi-process manner, but somehow cannot archieve the goal

def run(args):
    start = time.perf_counter()
    input_file = args.input
    chunk = args.chunk_size

    def cal_breaking(data):
        for index, row in data.iterrows():
            if row[1] == 0:  # mapping to the foward strand
                data.at[index, 'breaking_pos'] = int(row[5]) + int(row[3])
            elif row[1] == 16:  # mapping to the reverse strand
                data.at[index, 'breaking_pos'] = int(row[3])
            else:
                pass
        return data

    processes = []
    new_df = pd.DataFrame(columns=['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL'])
    for df in pd.read_csv(input_file, delimiter='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], chunksize=chunk):
        df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
        df = df.loc[~df['CIGAR'].str.contains('S') & ~df['CIGAR'].str.contains('H')]
        df['CIGAR'] = df.iloc[:, 5].str.extract(r'(\d+)')  # -d+ regex expression representing one or more numbers(0-9)
        df['breaking_pos'] = None
        with concurrent.futures.ProcessPoolExecutor(max_workers=6) as executor:
            processes.append(executor.submit(cal_breaking, df))

    for process in processes:
        new_df = pd.concat([new_df, process.result()], sort=True)

    new_df['count'] = 1
    new_df = new_df.groupby(['RNAME', 'breaking_pos']).count()['count'].reset_index()
    new_df['end'] = new_df['breaking_pos'] + 1
    new_df = new_df[['RNAME', 'breaking_pos', 'end', 'count']]
    print(new_df)
    end = time.perf_counter()
    print(f'process finished in {round(end - start, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    # parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.add_argument("-ck", help="read in chunk size", dest="chunk_size", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
