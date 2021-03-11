#! /usr/bin/env python

import argparse
import csv
import time
import gzip


# filter compressed pair file
def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    threshold = args.mq_threshold
    total_pair_count = 0
    qualified_pair_count = 0
    with gzip.open(input_file, 'rt') as input_pair:
        with gzip.open(output_file, 'wt') as output_pair:
            csv_reader = csv.reader(input_pair, delimiter='\t')
            csv_writer = csv.writer(output_pair, delimiter='\t')
            for line in csv_reader:
                if (int(line[8]) >= threshold) & (int(line[9]) >= threshold):
                    csv_writer.writerow(line)
                    qualified_pair_count += 1
                total_pair_count += 1
    end = time.perf_counter()

    print(f'{qualified_pair_count}  pairs has MPQ over {threshold}, taking {100*qualified_pair_count/total_pair_count}')
    print(f'process takes {round(end-start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="HiC pair MPQ filtering from compressed pair file")
    parser.add_argument("-i", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-t", help="MQ Threshold", dest="mq_threshold", type=int, required=True)
    parser.add_argument("-o", help="output pairs file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()