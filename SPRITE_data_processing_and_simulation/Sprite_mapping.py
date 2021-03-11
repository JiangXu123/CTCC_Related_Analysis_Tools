#! /usr/bin/env python

import pandas as pd
import argparse
import time
import pysam


def run(args):
    start = time.perf_counter()
    barcode_file = args.bc  # seq_file_R1 contains the DPM barcode and the genomic DNA sequencing information
    bam_file = args.bam  # the barcode-tagged reads number tsv file
    output_file = args.output
    df = pd.read_csv(barcode_file, names=['read_name', 'DPM_barcode', 'terminal_barcode', 'sec_odd_barcode', 'first_odd_barcode', 'even_barcode', 'complex_barcode'], delimiter='\t')
    bf = pysam.AlignmentFile(bam_file, 'rb')  # read in bam file indicated by 'rb'
    line = bf.fetch(until_eof=True)
    ls = []
    for i in range(0, len(df)):
        bam_line = next(line)
        read_name_1 = bam_line.qname
        read_name_2 = df.iloc[i, 0][1:].split(" ")[0]  # get the read name off the @ at the beginning, this name is also extracted by samtools
        while read_name_1 != read_name_2:  # some reads may have secondary mapping sites, or not mapped to the reference genome
            bam_line = next(line)
            read_name_1 = bam_line.qname
        map_q = bam_line.mapq
        read_chr = bam_line.reference_name
        if not bam_line.is_reverse:
            read_break_pos = bam_line.reference_start + 1
        if bam_line.is_reverse:
            read_break_pos = bam_line.reference_start + bam_line.query_alignment_length + 1
        if read_name_1 == read_name_2:  # \d+ means 1 or more number
            ls.append([read_name_1, read_chr, read_break_pos, map_q, df.iloc[i][6]])

    new_df = pd.DataFrame(ls, columns=['QNAME', 'mapped_chr', 'mapped_pos', 'MAPQ', 'barcode'])
    end = time.perf_counter()
    print(f'process {len(df)} reads takes {round(end - start, 2)} seconds')
    new_df.to_csv(output_file, header=False, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="Integrate read mapping position and barcode together")
    parser.add_argument("-rbc", help="read demultiplexed barcode file", dest="bc", type=str, required=True)
    parser.add_argument("-bam", help="mapped bam file", dest="bam", type=str, required=True)
    parser.add_argument("-out", help="the mapped barcode file ", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()