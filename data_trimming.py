#! /usr/bin/env python

import csv
import random as ran
import argparse
import time


def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    desired_length = args.des_len
    ori_pair_len = 0
    with open(input_file, 'r') as file1:
        for _ in file1:
            ori_pair_len += 1
    trim_len = ori_pair_len - desired_length
    print(f'{trim_len} pairs need to be randomly trimmed')
    rand_num_ls = ran.sample(range(0, ori_pair_len), trim_len+1)  # generate a random number list for reference
    rand_num_ls.sort()  # sort this random number list according to values from low to high
    print(f'{len(rand_num_ls)} non repeating random number generated')
    line_num = 1
    final_line_number = 0
    deleted_line_number = 0
    with open(input_file, 'r') as file1:
        with open(output_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            for line in csv_reader:
                if len(rand_num_ls) != 0:
                    if line_num != rand_num_ls[0]:
                        csv_writer.writerow(line)
                        final_line_number += 1
                    elif line_num == rand_num_ls[0]:
                        del rand_num_ls[0]
                        deleted_line_number += 1
                if len(rand_num_ls) == 0:
                    csv_writer.writerow(line)
                    final_line_number += 1
                line_num += 1
                if line_num % 10000 == 0:
                    print(f'{line_num} pairs processed')
    end = time.perf_counter()
    print(f'Original pair file randomly trimmed to {final_line_number} pairs in {round(end - start, 2)} seconds')
    print(f'{deleted_line_number} pairs are deleted')


def main():
    parser = argparse.ArgumentParser(description="randomly trimming the hic pairs file to the desired length")
    parser.add_argument("-i", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-o", help="output files name", dest="output", type=str, required=True)
    parser.add_argument("-l", help="desired final length", dest="des_len", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
