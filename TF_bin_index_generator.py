#! /usr/bin/env python

import argparse
import csv
import time
import os
import subprocess
# calculate the center locations of peaks from a batch of optimal IDR narrow peak files and make peak position file
# sort the peak position file using bash terminal sort -k1,1d k2,2n command
# make bins from the peak and use bgzip and tabix to index the bin file


def run(args):
    start = time.perf_counter()
    input_folder = args.input
    bin_radius = args.radius
    os.mkdir('./TF_peak_pos_output')
    input_file_name_ls = os.listdir(input_folder)
    input_file_path_ls = []
    peak_output_path_ls = []
    bin_folder_path = './TF_bin_' + str(bin_radius)
    os.mkdir(bin_folder_path)
    bin_file_path_ls = []
    compressed_bin_file_path_ls = []
    sorted_peak_file_path_ls = []

    i = 1
    TF_number_reference_ls = []
    for file_name in input_file_name_ls:
        input_file_path = input_folder + '/' + file_name
        input_file_path_ls.append(input_file_path)
        output_file_path = './TF_peak_pos_output' + '/' + file_name.split('_')[0] + '_' + 'peak_pos.tsv' #str(bin_radius) + 'bin'
        peak_output_path_ls.append(output_file_path)
        bin_file_path = bin_folder_path + '/' + str(i) + '_' + str(bin_radius) + 'bin'
        bin_file_path_ls.append(bin_file_path)
        compressed_bin_file = bin_file_path + '.gz'
        compressed_bin_file_path_ls.append(compressed_bin_file)
        TF_number_reference_ls.append([i, file_name.split('_')[0], file_name.split('.')[0]])
        i += 1

    with open('./TF_number_reference.tsv', 'w') as file1:  # store the TF name and
        csv_writer = csv.writer(file1, delimiter='\t')
        for item in TF_number_reference_ls:
            csv_writer.writerow(item)

    for i in range(0, len(input_file_path_ls)):
        input_file = input_file_path_ls[i]
        output_file = peak_output_path_ls[i]
        cal_peak_pos(input_file, output_file)

    for i in range(0, len(input_file_name_ls)):
        input_file = peak_output_path_ls[i]
        output_file = './TF_peak_pos_output/' + input_file_name_ls[i].split('_')[0] + '_peak_sorted.tsv'
        subprocess.call(f'sort -k1,1d -k2,2n {input_file} > {output_file}', shell=True)  # sort the peak file
        sorted_peak_file_path_ls.append(output_file)

    for i in range(0, len(input_file_name_ls)):
        input_file = sorted_peak_file_path_ls[i]
        output_file = bin_file_path_ls[i]
        bin_generator(input_file, output_file, bin_radius)
        compressed_bin_file = compressed_bin_file_path_ls[i]
        subprocess.call(f'sort -k1,1d -k2,2n {output_file} | bgzip > {compressed_bin_file}', shell=True)
        subprocess.call(f'tabix -s 1 -b 2 -e 3 {compressed_bin_file}', shell=True)

    end = time.perf_counter()
    print(f'process takes {round(end - start, 2)} seconds')


def cal_peak_pos(in_put, out_put):
    with open(in_put, 'r') as file1:
        with open(out_put, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            for line in csv_reader:
                peak_pos = int(line[1]) + int((int(line[2]) - int(line[1]))/2)
                csv_writer.writerow([line[0], peak_pos-1, peak_pos])


def bin_generator(in_put, out_put, bin_radius):
    pre_chr = 'chr1'
    with open(in_put, 'r') as peak_pos:
        with open(out_put, 'w') as bin_file:
            csv_writer = csv.writer(bin_file, delimiter='\t')
            csv_reader = csv.reader(peak_pos, delimiter='\t')
            pre_peak_pos = int(next(csv_reader)[2])
            i = 1
            for line in csv_reader:
                if pre_chr == line[0]:
                    if int(line[2]) >= pre_peak_pos + 2 * bin_radius:  # to guarantee there's no bin overlap
                        pre_peak_pos = int(line[2])
                        left_boarder = (int(line[2]) - bin_radius)
                        right_boarder = (int(line[2]) + bin_radius)
                        csv_writer.writerow([line[0], left_boarder, right_boarder, int(line[2]), str(i)])  # chr left_boarder right_boarder
                        i += 1  # chr_name, left boarder, right_boarder, center_pos,bin_number
                elif pre_chr != line[0]:  # if a new chr is encountered
                    pre_chr = line[0]
                    pre_peak_pos = int(line[2])


def main():
    parser = argparse.ArgumentParser(description="generate tabix indexed bin files from narrow peak bed file")
    parser.add_argument("-i", help="input narrow peak bed file folder", dest="input", type=str, required=True)
    parser.add_argument("-r", help="bin radius", dest="radius", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()