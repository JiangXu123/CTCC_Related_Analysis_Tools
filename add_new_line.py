#! /usr/bin/env python

import argparse
import csv
import time


def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    new_line_contents = args.new_l_c

    with open(input_file, 'r') as file1:
        with open(output_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            csv_writer.writerow([new_line_contents, 0, 0, 0])
            for line in csv_reader:
                csv_writer.writerow(line)
    end = time.perf_counter()
    print(f'process takes {round(end-start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="add a new line to the beginning of a tsv file")
    parser.add_argument("-i", help="input file", dest="input", type=str, required=True)
    parser.add_argument("-o", help="output file", dest="output", type=str, required=True)
    parser.add_argument("-c", help="contents of the new line", dest="new_l_c", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()