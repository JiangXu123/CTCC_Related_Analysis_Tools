#! /usr/bin/env python

import csv
import gzip
import argparse
import time


def run(args):
    changed_name = args.rname
    input_file = args.input
    output_file = args.output
    i = 1
    start = time.perf_counter()
    with gzip.open(input_file, 'rt') as file1:  # open and write in gzip compressed text format
        with gzip.open(output_file, 'wt') as file2:
            csv_reader = csv.reader(file1, delimiter=' ')
            csv_writer = csv.writer(file2, delimiter=' ')
            for line in csv_reader:
                if i % 4 == 1:
                    name_1 = changed_name + '.' + str(round(i / 4) + 1)
                    name_2 = line[1]
                    csv_writer.writerow([name_1, name_2])  # each read name is separated into two part
                else:
                    csv_writer.writerow(line)
                i += 1
    end = time.perf_counter()
    print(f'{int((i-1)/4)}  reads name corrected, taking {round(end - start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="repair the read name for correct HiC-Pro processing")
    parser.add_argument("-i", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-o", help="output pairs file", dest="output", type=str, required=True)
    parser.add_argument("-rn", help="replaced name", dest="rname", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
