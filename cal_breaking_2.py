#! /usr/bin/env python

import pandas as pd
import argparse
import time
import tabix

def run(args):
    start1 = time.perf_counter()
    input_file = args.input
    chunk = args.chunk_size
    output_file = args.output
    trim_sel = args.trim_s
    mq = args.mpq
    new_df = pd.DataFrame(columns=['RNAME', 'breaking_pos', 'count'])
    idx = 0
    for df in pd.read_csv(input_file, delimiter='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], chunksize=chunk, header=None):
        start = time.perf_counter()
        df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ', 'QUAL']
        df = df.loc[(~df['CIGAR'].str.contains('S')) & (~df['CIGAR'].str.contains('H')) & (~df['CIGAR'].str.contains('I')) & (~df['CIGAR'].str.contains('D')) & (~df['CIGAR'].str.contains('P')) & (~df['CIGAR'].str.contains('N')) & (~df['CIGAR'].str.contains('X')) & (df['MAPQ'] >= mq)]
        df['CIGAR'] = df.iloc[:, 5].str.extract(r'(\d+)')
        for index, row in df.iterrows():
            if trim_sel == 3:  # connector is trimmed away from 3'
                if row[1] == 0:  # mapping to the forward strand
                    df.at[index, 'breaking_pos'] = (int(row[5]) + int(row[3]) - 1)
                elif row[1] == 16:  # mapping to the reverse strand
                    df.at[index, 'breaking_pos'] = int(row[3])
                else:
                    pass
            if trim_sel == 5:  # connector is trimmed away from 3'
                if row[1] == 0:  # mapping to the forward strand
                    df.at[index, 'breaking_pos'] = int(row[3])
                elif row[1] == 16:  # mapping to the reverse strand
                    df.at[index, 'breaking_pos'] = (int(row[5]) + int(row[3]) - 1)
                else:
                    pass
        df['count'] = 1
        df = df[['RNAME', 'breaking_pos', 'count']]
        end = time.perf_counter()
        new_df = pd.concat([new_df, df], sort=True)
        print(f'finish the {idx} {chunk} chunk in {round(end - start, 2)} sec(s)')
        idx += 1
    new_df = new_df.groupby(['RNAME', 'breaking_pos']).count()['count'].reset_index()
    new_df['start'] = new_df['breaking_pos'] - 1
    new_df = new_df[['RNAME', 'start', 'breaking_pos', 'count']]
    new_df[['start', 'breaking_pos']] = new_df[['start', 'breaking_pos']].astype(int)
    new_df.to_csv(output_file, '\t', index=None, header=None)
    end1 = time.perf_counter()
    print(f'process finished in {round(end1 - start1, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="compute breaking point of CTCC or HiC data based on adaptor trimmed "
                                                 "mapping sam file")
    parser.add_argument("-in", help="input sam file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output tsv bed file", dest="output", type=str, required=True)
    parser.add_argument("-mqt", help="MAPQ threshold", dest="mpq", type=int, required=True)
    parser.add_argument("-ck", help="read in chunk size", dest="chunk_size", type=int, required=True)
    parser.add_argument("-trim", help="connector trimming fassion 3 or 5 ", dest="trim_s", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
