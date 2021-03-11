#! /usr/bin/env python

import argparse
import tabix
import csv
import time


def run(args):
    start = time.perf_counter()
    seg_start = time.perf_counter()
    compartment_filename = args.comp
    input_pair_file = args.input
    output_pair_compartment_file = args.output
    experiment = args.exp
    tb = tabix.open(compartment_filename)
    line_count = 0
    valid_line_count = 0
    with open(input_pair_file, 'r') as pairs_file:
        with open(output_pair_compartment_file, 'w') as pairs_compartment_file:
            csv_reader = csv.reader(pairs_file, delimiter='\t')
            csv_writer = csv.writer(pairs_compartment_file, delimiter='\t')
            for line in csv_reader:
                try:
                    result1 = tb.query(line[1], int(line[2]) - 1, int(line[2]))
                    result2 = tb.query(line[4], int(line[5]) - 1, int(line[5]))
                    read1_comp = next(result1)[3]
                    read2_comp = next(result2)[3]
                    interaction_type = read1_comp + '-' + read2_comp
                    csv_writer.writerow([line[1], line[2], line[3], line[4], line[5], line[6], interaction_type, experiment])
                    valid_line_count += 1
                except:
                    pass
                line_count += 1
                if line_count%1000000 == 0:
                    seg_end = time.perf_counter()
                    print(f'{line_count} pairs processed, takes {seg_end - seg_start} sec')
                    seg_start = 0
    end = time.perf_counter()
    print(f'from {line_count} pairs, {valid_line_count} pairs tagged in {round(end-start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.add_argument("-exp", help="experiment", dest="exp", type=str, required=True)
    parser.add_argument("-cf", help="compartment file", dest="comp", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
