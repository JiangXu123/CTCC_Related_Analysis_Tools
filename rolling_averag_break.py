#! /usr/bin/env python
# this python code is used to extract the bed file according to the region in the genome: RNAME, start, end.
import pandas as pd
import time
import argparse
import tabix


def run(args):
    start1 = time.perf_counter()
    window = args.window_size
    chr_sel = args.chrsel
    chr_start = args.start
    chr_end = args.end
    slide_step = args.slidestep
    rf_bf = args.reference
    chromsize_file = args.chrsf
    chrom_size = pd.read_csv(chromsize_file, delimiter='\t', header=None)
    for index, row in chrom_size.iterrows():
        if row[0] == chr_sel:                  # read the size of the selected the chromosome,
            size = row[1]                      # which is contained in the second column
    # if ~chr_end:           # if no chr_end specified
    #     chr_end = size
    # if ~chr_start:
    #     chr_start = 1	   # if no chr_start specified
    tb = tabix.open(rf_bf)
    new_ls = []
    for i in range(chr_start, chr_end - window + 1, slide_step):
        counts = 0
        tb_results = tb.query(chr_sel, i, i+window-1)
        for result in tb_results:
            counts += int(result[3])   # result[3] is the break counts in that particular position
        new_ls.append([chr_sel, i+round(window/2)-1, i+round(window/2), counts])  # row[3] is the position of the center of the window
        print(i)
    new_df = pd.DataFrame(new_ls, columns=['RNAME', 'Center_pos', 'end', 'average_counts'])
    new_df.to_csv(f'{chr_sel}_{chr_start}_{chr_end}_averaged_{window}_breaking_counts.bed', header=False, index=False, sep='\t')
    end1 = time.perf_counter()
    print(f'process finished in {round(end1 - start1, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="generating a sliding window average bedfile for breaking data")
    parser.add_argument("-rb", help="Reference bed file that contains the breaking information, compressed bgzip and indexed with tabix", dest="reference", type=str, required=True)
    parser.add_argument("-chrs", help="chromesome size file", dest="chrsf", type=str, required=True)
    parser.add_argument("-chr_sel", help="select the chromosome", dest="chrsel", type=str, required=True)
    parser.add_argument("-winds", help="sliding window size", dest="window_size", type=int, required=True)
    parser.add_argument("-sel_stt", help="select range, start", dest="start", type=int, required=False)
    parser.add_argument("-sel_end", help="select range, end", dest="end", type=int, required=False)
    parser.add_argument("-stp", help="sliding step", dest="slidestep", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()