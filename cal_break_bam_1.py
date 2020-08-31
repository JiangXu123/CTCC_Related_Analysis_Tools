#! /usr/bin/env python
# This program is used to calculate aggregated breaking point information from MNase_seq or SPRITE experiment
import pysam
import time
import argparse
import pandas as pd


def run(args):
    start1 = time.perf_counter()
    bam_file = args.input
    output_file = args.output
    mq = args.mpq
    bf = pysam.AlignmentFile(bam_file, 'rb')  # use the 'until_eof=True' option to iterate through the file
    count = 1
    new_df = pd.DataFrame(columns=['RNAME', 'start', 'count'])
    i = 0
    s = 0
    h = 0
    d = 0
    p = 0
    n = 0
    x = 0
    w = 0   # total mapping counts that equal or exceed the mapq threshold
    tw = 0  # total whole mappings
    for line in bf.fetch(until_eof=True):
        if line.mapq >= mq:
            if 'I' in line.cigarstring:
                i += 1
            if 'S' in line.cigarstring:
                s += 1
            if 'H' in line.cigarstring:
                h += 1
            if 'D' in line.cigarstring:
                d += 1
            if 'P' in line.cigarstring:
                p += 1
            if 'N' in line.cigarstring:
                n += 1
            if 'X' in line.cigarstring:
                x += 1
            else:
                if line.flag == 0:
                    df = pd.DataFrame([[line.reference_name, line.reference_start, count]],
                                  columns=['RNAME', 'start', 'count'])
                if line.flag == 16:
                    df = pd.DataFrame([[line.reference_name, line.reference_start + line.cigar[0][1], count]],
                                  columns=['RNAME', 'start', 'count'])
                new_df = pd.concat([new_df, df], sort=False)
            W += 1
        TW += 1
    new_df.reset_index(inplace=True, drop=True)
    new_df = new_df.groupby(['RNAME', 'start']).count()['count'].reset_index()
    new_df['breaking_pos'] = new_df['start'] + 1
    new_df = new_df[['RNAME', 'start', 'breaking_pos', 'count']]
    new_df.to_csv(output_file, sep='\t', header=False, index=False)
    end1 = time.perf_counter()
    print(f'mapping with insertion(I) take {100 * I / W} percent')
    print(f'mapping with insertion(I) has {I} count')
    print(f'mapping with deletion(D) take {100 * D / W} percent')
    print(f'mapping with deletion(D) take {D} count')
    print(f'mapping with soft clipping(S) take {100 * S / W} percent')
    print(f'mapping with soft clipping(S) take {S} count')
    print(f'mapping with hard clipping(H) take {100 * H / W} percent')
    print(f'mapping with hard clipping(H) take {H} count')
    print(f'mapping with padding(P) take {100 * P / W} percent')
    print(f'mapping with padding(P) take {100 * P / W} percent')
    print(f'mapping with skipped bases(N) take {100 * N / W} percent')
    print(f'mapping with mismatch(X) take {100 * X / W} percent')
    print(f'process finished in {round(end1 - start1, 2)} second(s)')
    print(f'mappings with mapq equal or more than {mq} take {100 * W/TW} percent')


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-in", help="input bamfile", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output bed file name", dest="output", type=str, required=True)
    parser.add_argument("-mqt", help="MAPQ threshold", dest="mpq", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()