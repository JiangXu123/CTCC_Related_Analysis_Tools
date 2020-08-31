#! /usr/bin/env python

import argparse
import time
import os

# concatenate the files in a folder


def run(args):
    start = time.perf_counter()
    input_folder = args.infolder  # folder path containing the files to be concatenated
    output_file = args.otf
    input_file_name_ls = os.listdir(input_folder)
    input_file_path_ls = []

    for file in input_file_name_ls:
        input_file_path = os.path.abspath(input_folder) + '/' + file
        input_file_path_ls.append(input_file_path)

    with open(output_file, 'w') as file1:
        for input_file in input_file_path_ls:
            with open(input_file, 'r') as file2:
                for line in file2:
                    file1.write(line)
    end = time.perf_counter()
    print(f'{len(input_file_name_ls)} file(s) processed within {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="concatenate the files in a folder")
    parser.add_argument("-i", help="peak signal centered indexed bin file", dest="infolder", type=str, required=True)
    parser.add_argument("-o", help="output file name", dest="otf", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
