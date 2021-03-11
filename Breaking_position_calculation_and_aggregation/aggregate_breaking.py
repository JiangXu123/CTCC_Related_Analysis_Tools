#! /usr/bin/env python

import csv
import time
import argparse


def run(args):
    start1 = time.perf_counter()
    input_file = args.input
    output_file = args.output
    total_breaks = 0
    with open(input_file, 'r') as file1:
        for _ in file1:
            total_breaks += 1
    print(f"total {total_breaks} breaking point found")
    aggr_count = 1
    with open(input_file, 'r') as file1:
        with open(output_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            temp_line = next(csv_reader)
            temp_break_pos = temp_line[2]
            temp_chr = temp_line[0]
            for line in csv_reader:
                if (temp_break_pos == line[2]) & (temp_chr == line[0]):
                    aggr_count += 1
                elif (temp_break_pos != line[2]) | (temp_chr != line[0]):
                    csv_writer.writerow([line[0], int(line[1])-1, int(line[2])-1, aggr_count, round(1000000*aggr_count/total_breaks, 5)])
                    temp_break_pos = line[2]
                    aggr_count = 1
    end1 = time.perf_counter()
    print(f'process finished in {round(end1 - start1, 2)} second(s)')


def main():
    parser = argparse.ArgumentParser(description="aggregate the breaking counts to breakings per million breakings")
    parser.add_argument("-i", help="input breaking position file with compartment information, should be sorted with sort -k1,1a, -k2,2n", dest="input", type=str, required=True)
    parser.add_argument("-o", help="output aggregated file name", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()