#! /usr/bin/env python

import csv
import argparse
import time
import tabix
import pysam


# to extract breaking position information and subcompartment information from a CTCC sequencing result
def run(args):
    start = time.perf_counter()
    compartment_filename = args.comp
    bam_file = args.input
    breaking_compartment_tagged_file = args.output
    trim_sel = args.trim_s
    mq = args.mpq
    chromosome_ls = []
    for i in range(1, 23):  # generate a list that contains all the chromosomes' names that the compartment file has
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')
    chromosome_ls.append('chrY')
    chromosome_ls.append('chrM')
    tb = tabix.open(compartment_filename)
    line_count = 0
    valid_count = 0
    print(f'chromosomes searched are {chromosome_ls}')
    bf = pysam.AlignmentFile(bam_file, 'rb')
    with open(breaking_compartment_tagged_file, 'w') as breaking_compartment_file:
        csv_writer = csv.writer(breaking_compartment_file, delimiter='\t')
        for line in bf.fetch(until_eof=True):
            if line.mapq >= mq:
                if trim_sel == 3:  # connector is trimmed away from 3'
                    if not line.is_reverse:  # if mapped to the forward strand
                        breaking_pos = line.reference_start + line.query_alignment_length - 1
                        chr_name = line.reference_name
                    if line.is_reverse: # if mapped to the reverse strand
                        breaking_pos = line.reference_start
                        chr_name = line.reference_name
                if trim_sel == 5:  # connector is trimmed away from 3'
                    if line.is_reverse:  # if mapped to the reverse strand
                        breaking_pos = line.reference_start + line.query_alignment_length - 1
                        chr_name = line.reference_name
                    if not line.is_reverse:  # if mapped to the forward strand
                        breaking_pos = line.reference_start
                        chr_name = line.reference_name
                if chr_name in chromosome_ls: # chromosome_ls contains all the chromosomes' names that the compartment file has
                    result1 = tb.query(chr_name, int(breaking_pos) - 1, int(breaking_pos)) # get the subcompartment information
                    try:  # if the breaking position is contained in the subcompartment file
                        breaking_pos_comp = next(result1)[3]
                    except: # if the breaking position is not contained in the subcompartment file
                        breaking_pos_comp = None
                if breaking_pos_comp:   # if it return a none 'None' value
                    csv_writer.writerow([chr_name, breaking_pos-1, breaking_pos, breaking_pos_comp])
                valid_count += 1   # valid_count contains number of the breaking position that has mpq >= mq
            line_count += 1
    end = time.perf_counter()
    print(f'from {line_count} total breaking point tagged {valid_count} breaking point with mapping quality >= {mq} in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="extract breaking position and subcompartment information from adaptor trimmed CTCC bam file")
    parser.add_argument("-in", help="input sam file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output tsv bed file", dest="output", type=str, required=True)
    parser.add_argument("-mqt", help="MAPQ threshold", dest="mpq", type=int, required=True)
    parser.add_argument("-cp", help="tabix indexed bgziped compartment file", dest="comp", type=str, required=True)
    parser.add_argument("-trim", help="connector trimming fashion 3 or 5 ", dest="trim_s", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
