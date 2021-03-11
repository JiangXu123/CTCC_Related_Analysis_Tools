#! /usr/bin/env python

import argparse
import csv
import os
import shutil
# rename files to according numerical sequence and make a tsv file contain the old file
# and new file name


def run(args):
    input_folder = args.input
    rename_file_folder = input_folder + '_renamed'
    os.mkdir(rename_file_folder)
    old_file_name_ls = os.listdir(input_folder)
    old_file_path_ls = []
    for file in old_file_name_ls:
        old_file_path = input_folder + '/' + file
        old_file_path_ls.append(old_file_path)
    old_file_suffix = old_file_name_ls[0].split('.')[-1]
    new_file_name_index = 0
    file_name_reference = input_folder + '_rename-reference.tsv'
    with open(file_name_reference, 'w') as file1:
        csv_writer = csv.writer(file1, delimiter='\t')
        for file in old_file_path_ls:
            old_file_size = os.path.getsize(file)
            new_file_path = rename_file_folder + '/' + f'{new_file_name_index}.'+old_file_suffix
            shutil.copy(file, new_file_path)
            new_file_size = os.path.getsize(new_file_path)
            if old_file_size == new_file_size:
                print(f'renamed file {new_file_name_index} size checked correct')
                new_file_name = f'{new_file_name_index}.'+old_file_suffix
                csv_writer.writerow([file.split('/')[-1], new_file_name, new_file_size])  # write the reference list
                new_file_name_index += 1


def main():
    parser = argparse.ArgumentParser(description="rename file in numerical order to make it easier to process in batch")
    parser.add_argument("-i", help="folder that contain the file to be renamed", dest="input", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
