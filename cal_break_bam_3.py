#! /usr/bin/env python
# This program is used to calculate aggregated breaking point information from MNase_seq or SPRITE experiment
import pysam
import time
import argparse
import csv


def run(args):
    start1 = time.perf_counter()
    bam_file = args.input
    output_file = args.output
    mq = args.mpq
    bf = pysam.AlignmentFile(bam_file, 'rb')  # use the 'until_eof=True' option to iterate through the file
    breaking_posi_ls = []
    i = 0
    s = 0
    h = 0
    d = 0
    p = 0
    n = 0
    x = 0
    w = 0   # total mapping counts that equal or exceed the mapq threshold
    tw = 0  # total whole mappings
    with open(output_file, 'w') as file:
        csv_writer = csv.writer()
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
                if not line.is_reverse:   # if mapping to the forward strand
                    breaking_posi_ls.append([line.reference_name, line.reference_start + 1, 1])
                if line.is_reverse:
                    breaking_posi_ls.append([line.reference_name, line.reference_start + line.query_alignment_length + 1, 1])
                w += 1
            tw += 1
    new_df = pd.DataFrame(breaking_posi_ls, columns=['RNAME', 'start', 'count'])
    new_df = new_df.groupby(['RNAME', 'start']).count()['count'].reset_index()
    new_df['breaking_pos'] = new_df['start'] - 1
    new_df = new_df[['RNAME', 'breaking_pos', 'start', 'count']]
    new_df.to_csv(output_file, sep='\t', header=False, index=False)
    end1 = time.perf_counter()
    print(f'mapping with insertion(I) take {100 * i / w} percent')
    print(f'mapping with insertion(I) has {i} count')
    print(f'mapping with deletion(D) take {100 * d / w} percent')
    print(f'mapping with deletion(D) take {d} count')
    print(f'mapping with soft clipping(S) take {100 * s / w} percent')
    print(f'mapping with soft clipping(S) take {s} count')
    print(f'mapping with hard clipping(H) take {100 * h / w} percent')
    print(f'mapping with hard clipping(H) take {h} count')
    print(f'mapping with mismatch(X) take {100 * x / w} percent')
    print(f'process finished in {round(end1 - start1, 2)} second(s)')
    print(f'mappings with mapq equal or more than {mq} take {100 * w/tw} percent')


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