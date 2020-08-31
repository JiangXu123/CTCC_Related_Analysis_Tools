#! /usr/bin/env python

import argparse
import tabix
import csv
import time


def run(args):
    start = time.perf_counter()
    compartment_filename = args.comp
    breaking_pos_file = args.input
    breaking_compartment_tagged_file = args.output
    chromosome_ls = []
    for i in range(1, 23):
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')
    chromosome_ls.append('chrY')
    chromosome_ls.append('chrM')
    tb = tabix.open(compartment_filename)
    line_count = 0
    print(f'chromosomes searched are {chromosome_ls}')
    with open(breaking_pos_file, 'r') as breaking_file:
        with open(breaking_compartment_tagged_file, 'w') as breaking_compartment_file:
            csv_reader = csv.reader(breaking_file, delimiter='\t')
            csv_writer = csv.writer(breaking_compartment_file, delimiter='\t')
            for line in csv_reader:
                if line[0] in chromosome_ls:
                    result1 = tb.query(line[0], int(line[1]), int(line[2]))
                    breaking_pos_comp = next(result1)[3]
                    csv_writer.writerow([line[0], int(line[1]), int(line[2]), int(line[3]), breaking_pos_comp])
                    line_count += 1
    end = time.perf_counter()
    print(f'{line_count} pairs tagged in {round(end-start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="subcompartment tagging of breaking position")
    parser.add_argument("-in", help="input breaking position file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output subcompartment tagged file", dest="output", type=str, required=True)
    parser.add_argument("-cf", help="compartment file", dest="comp", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
