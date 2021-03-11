#! /usr/bin/env python

import argparse
import gzip
import time


# trim connector sequence from 5'
def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    trimming_number = args.num_trim
    line_number = 1
    with gzip.open(input_file, 'rt') as input_file_obj:
        with gzip.open(output_file, 'wt') as output_file_obj:
            for line in input_file_obj:
                if (line_number % 4 == 2) | (line_number % 4 == 0):  # both sequence line and quality line should be trimmed
                    output_file_obj.write(line[trimming_number:])
                else:
                    output_file_obj.write(line)
                line_number += 1
    end = time.perf_counter()
    print(f'{trimming_number} bases trimmed from 5 prime end for {line_number/4} reads')
    print(f'process takes {round(end-start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="trimming reads from 5' end by specifying the trimming length")
    parser.add_argument("-i", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-n", help="number of bases trimmed", dest="num_trim", type=int, required=True)
    parser.add_argument("-o", help="output pairs file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()