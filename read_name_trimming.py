#! /usr/bin/env python

import argparse
import csv
import time
import gzip


def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    k = 1
    with gzip.open(input_file, 'rt') as file1:
        with gzip.open(output_file, 'wt') as file2:
            csv_reader = csv.reader(file1, delimiter=' ')
            csv_writer = csv.writer(file2, delimiter=' ')
            for line in csv_reader:
                if k % 4 == 1:
                    temp_ls = line[0].split('.')[0:3]
                    temp_name = ['.'.join(temp_ls[0:3])]
                    csv_writer.writerow(temp_name)
                if k % 4 == 2:
                    csv_writer.writerow(line)
                if k % 4 == 3:
                    temp_ls = line[0].split('.')[0:3]
                    temp_name = ['.'.join(temp_ls[0:3])]
                    csv_writer.writerow(temp_name)
                if k % 4 == 0:
                    csv_writer.writerow(line)
                k += 1
    end = time.perf_counter()
    print(f'{k/4}  reads name corrected, taking {round(end-start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="HiC pair MPQ filtering from compressed pair file")
    parser.add_argument("-i", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-o", help="output pairs file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()