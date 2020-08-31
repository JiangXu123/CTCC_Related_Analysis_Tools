#! /usr/bin/env python
# this python code is used to extract the bed file according to the region in the genome: RNAME, start, end.
import pandas as pd
import time
import argparse
import tabix


def run(args):
    start1 = time.perf_counter()
    query_bed = args.query_b
    chr_start = args.start
    chr_end = args.end
    rf_bf = args.reference
    tb = tabix.open(rf_bf)
    new_ls = []
    df = pd.read_csv(query_bed, delimiter='\t', names=['RNAME', 'start', 'end', 'central'])
    df = df.loc[(df['start'] >= chr_start) & (df['end'] <= chr_end)]
    for index, row in df.iterrows():
        tb_results = tb.query(row[0], row[1], row[2])
        for result in tb_results:
            counts += int(result[3])   # result[3] is the break counts in that particular position
        new_ls.append([row[0], row[3], row[3]+1, counts])  # row[3] is the position of the center of the window
    new_df = pd.DataFrame(new_ls, columns=['RNAME', 'Center_pos', 'end', 'average_counts'])
    new_df.to_csv(f'{chr_sel}_{chr_start}_{chr_end}_averaged_{window}_breaking_counts.bed', header=False, index=False, sep='\t')
    end1 = time.perf_counter()
    print(f'process finished in {round(end1 - start1, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="generating a sliding window average bedfile for breaking data")
    parser.add_argument("-rb", help="Reference bed file that contains the breaking information, compressed bgzip and indexed with tabix", dest="reference", type=str, required=True)
    parser.add_argument("-query_bed", help="select the chromosome", dest="query_b", type=str, required=True)
    parser.add_argument("-sel_stt", help="select range, start", dest="start", type=int, required=True)
    parser.add_argument("-sel_end", help="select range, end", dest="end", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()