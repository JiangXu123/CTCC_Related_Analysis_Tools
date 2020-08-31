#! /usr/bin/env python

import argparse
import pysam
import tabix
import csv
import time

# extract breaking information from the bam file of SPRITE, ChIP_Seq, MNase_Seq etc... that don't form ligation product
# and tag subcompartment information or any tabix indexed property
def run(args):
    start = time.perf_counter()
    compartment_filename = args.comp
    bam_file = args.input
    breaking_compartment_tagged_file = args.output
    chromosome_ls = []
    for i in range(1, 23):  # generate a list that contains all the chromosomes' names that the compartment file has
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')
    chromosome_ls.append('chrY')
    chromosome_ls.append('chrM')
    tb = tabix.open(compartment_filename)
    line_count = 0
    print(f'chromosomes searched are {chromosome_ls}')
    bf = pysam.AlignmentFile(bam_file, 'rb')
    with open(breaking_compartment_tagged_file, 'w') as breaking_compartment_file:
        csv_writer = csv.writer(breaking_compartment_file, delimiter='\t')
        for line in bf.fetch(until_eof=True):
            if not line.is_reverse:  # if mapped to the forward strand
                chr_name = line.reference_name
                breaking_pos = line.reference_start
            if line.is_reverse:  # if mapped to the reverse strand
                chr_name = line.reference_name
                breaking_pos = line.reference_start + line.query_alignment_length - 1
            if chr_name in chromosome_ls:
                result1 = tb.query(chr_name, int(breaking_pos) - 1, int(breaking_pos))
                try:
                    breaking_pos_comp = next(result1)[3]
                except:
                    breaking_pos_comp = None
            if breaking_pos_comp:           # if it return a none 'None' value
                csv_writer.writerow([chr_name, breaking_pos-1, breaking_pos, breaking_pos_comp])
            line_count += 1

    end = time.perf_counter()
    print(f'{line_count} pairs tagged in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="extract breaking information from the bam file of SPRITE, ChIP_Seq, MNase_Seq etc, that don't form ligation product, tag subcompartment information or any tabix indexed property")
    parser.add_argument("-in", help="input filtered bam file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output subcompartment tagged file", dest="output", type=str, required=True)
    parser.add_argument("-cf", help="compartment file", dest="comp", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
