#! /usr/bin/env python

import argparse
import csv
import time
import tabix
import shutil
import subprocess
import multiprocessing
import os

# extract breaking information from the breaking position with compartment file
# and calculate the distance to the center of the nucleosomes or a TF binding sites(location)
# aggregate the breaking distance to nucleosome center or TF binding center


def run(args):
    start = time.perf_counter()
    thread_number = args.t_number
    indexed_bin_file = args.ib
    bin_radius = args.radius
    input_folder = args.bpfd  # folder path containing the breaking position file
    temp_path = os.path.abspath('..') + '/temp'
    input_file_name_ls = os.listdir(input_folder)
    input_file_path_ls = []
    output_file_path_ls = []
    aggr_output_file_path_ls = []
    bin_name = indexed_bin_file.split('/')[-1].split('.')[0]
    for file in input_file_name_ls:
        input_file_path = os.path.abspath(input_folder) + '/' + file
        input_file_path_ls.append(input_file_path)
        output_file_path = './' + file.split('.')[0] + '_' + bin_name + '_dist.tsv'
        # the output_file_path is the path of intermediate output_file contain breaking_dist_to_ref information
        output_file_path_ls.append(output_file_path)
        aggr_output_file_path = './' + file.split('.')[0] + '_' + bin_name + '_dist_aggr.tsv'
        aggr_output_file_path_ls.append(aggr_output_file_path)
    for i in range(0, len(input_file_path_ls)):
        start_1 = time.perf_counter()
        input_file = input_file_path_ls[i]
        input_file_name = input_file_name_ls[i]
        output_file = output_file_path_ls[i]
        aggr_output_file = aggr_output_file_path_ls[i]
        try:
            os.makedirs(temp_path)  # make a temporary folder in the current directory
            print(f'The temporary folder created at {temp_path}')
        except:
            print('temp folder already exist,existing file will be overwritten.')
        valid_line_count = 0
        total_line_count = 0
        with open(input_file, 'r') as file1:
            for _ in file1:
                total_line_count += 1  # count the total number of breakings

        experiment = input_file_name.split('-')[0]  # file name should be writen in 'experiment-**.tsv' format

        block_lines = int(total_line_count / thread_number)  # make an estimation how many lines should contain per block
        print(f'total of {total_line_count} lines will be split into {block_lines} per temporary block file')

        os.chdir(temp_path)  # enter the temporary folder
        subprocess.call(f'cat {input_file} | split -d -l {block_lines} - temp_in_', shell=True)  # I cheated here using bash shell's split method, which is much faster in splitting files
        os.chdir('../..')  # get back from the temporary folder
        block_num = len(os.listdir('./temp')) # block_num contains the actual split file number
        print(f'Input file {i} is split to {block_num} blocks into the temp folder')
        # the above code create the temporary folder and splits the file into smaller files
        processes = []
        temp_in_file_name_ls = os.listdir('./temp')
        temp_in_file_path_ls = []
        temp_out_file_path_ls = []
        for name in temp_in_file_name_ls:
            temp_in_file_path_ls.append('./temp/' + name)
            temp_out_file_path_ls.append('./temp/temp_out_' + name.split('_')[-1])

        for temp_index in range(0, block_num):  # the zip method allow to open two list in parallel
            in_file = temp_in_file_path_ls[temp_index]
            out_file = temp_out_file_path_ls[temp_index]
            proc = multiprocessing.Process(target=process_file, args=[in_file, out_file, bin_radius, indexed_bin_file, experiment])
            proc.start()
            processes.append(proc)

        for process in processes:
            process.join()  # wait until all processes finished

        with open(output_file, 'w') as file1:  # concatenate the temporary output files
            for file in temp_out_file_path_ls:
                with open(file, 'r') as file2:
                    for line in file2:
                        file1.write(line)
                        valid_line_count += 1
        shutil.rmtree(temp_path)  # remove the temporary folder after finishing processing
        aggregate_nc_dist(output_file, aggr_output_file, experiment, bin_name, bin_radius)
        os.unlink(output_file)  # remove the intermediate 'breaking_dist_to_ref' file
        end_1 = time.perf_counter()
        print(f'temporary folder removed')
        print(f'for file {input_file_name_ls[i]}, from {total_line_count} breaking positions, {valid_line_count} relative positions calculated in  {round(end_1 - start_1, 2)} sec')
    end = time.perf_counter()
    print(f'{len(input_file_name_ls)} file(s) processed within {round(end - start, 2)} sec')


def process_file(in_file, out_file, bin_radius, indexed_peak_file, exp):
    with open(in_file, 'r') as break_pos:
        with open(out_file, 'w') as aggregated_data:
            csv_writer = csv.writer(aggregated_data, delimiter='\t')
            csv_reader = csv.reader(break_pos, delimiter='\t')
            tb = tabix.open(indexed_peak_file)
            for line in csv_reader:
                try:
                    result = tb.query(line[0], int(line[1]), int(line[2])) # result1 contains: chr_name, left boarder, right_boarder, height, center_pos,bin_number
                    item = next(result)
                    distance = int(line[2]) - int(item[4])  # item[4] contains the bin center position
                    if abs(distance) <= bin_radius:
                        csv_writer.writerow([exp, int(item[5]), line[3], line[0], distance])
                        # columns name are: experiment name, nucleosome number(bin number), sub_nucleus compartment, chr_name,
                        # distance to nucleosome center
                except:
                    pass


def aggregate_nc_dist(input_dist_file, aggr_output_file, experiment, ref_bins_name, bin_radius):
    start_2 = time.perf_counter()
    A1_dic = {}
    A2_dic = {}
    B1_dic = {}
    B2_dic = {}
    B3_dic = {}
    B4_dic = {}
    A1_bin_ls = []
    A2_bin_ls = []
    B1_bin_ls = []
    B2_bin_ls = []
    B3_bin_ls = []
    B4_bin_ls = []

    compartment_ls = ['A1', 'A2', 'B1', 'B2', 'B3', 'B4']
    line_number = 0

    for i in range(-bin_radius, bin_radius+1):
        A1_dic.update({i: 0})
        A2_dic.update({i: 0})
        B1_dic.update({i: 0})
        B2_dic.update({i: 0})
        B3_dic.update({i: 0})
        B4_dic.update({i: 0})

    with open(input_dist_file, 'r') as break_pos:
        with open(aggr_output_file, 'w') as aggregated_data:
            csv_reader = csv.reader(break_pos, delimiter='\t')
            csv_writer = csv.writer(aggregated_data, delimiter='\t')
            for line in csv_reader:
                try:
                    if line[2] == 'A1':
                        A1_dic[int(line[4])] += 1  # line[4] records the distance from nucleosome center.
                        A1_bin_ls.append(line[1])  # line[1] records the nucleosome ID number(each nucleosome has an ID number)
                    elif line[2] == 'A2':
                        A2_dic[int(line[4])] += 1
                        A2_bin_ls.append(line[1])
                    elif line[2] == 'B1':
                        B1_dic[int(line[4])] += 1
                        B1_bin_ls.append(line[1])
                    elif line[2] == 'B2':
                        B2_dic[int(line[4])] += 1
                        B2_bin_ls.append(line[1])
                    elif line[2] == 'B3':
                        B3_dic[int(line[4])] += 1
                        B3_bin_ls.append(line[1])
                    elif line[2] == 'B4':
                        B4_dic[int(line[4])] += 1
                        B4_bin_ls.append(line[1])
                except:
                    pass
                line_number += 1
            A1_bin_count = len(list(dict.fromkeys(A1_bin_ls)))
            A2_bin_count = len(list(dict.fromkeys(A2_bin_ls)))
            B1_bin_count = len(list(dict.fromkeys(B1_bin_ls)))
            B2_bin_count = len(list(dict.fromkeys(B2_bin_ls)))
            B3_bin_count = len(list(dict.fromkeys(B3_bin_ls)))
            B4_bin_count = len(list(dict.fromkeys(B4_bin_ls)))
            print(f'bin belong to A1 have {A1_bin_count}')
            print(f'bin belong to A2 have {A2_bin_count}')
            print(f'bin belong to B1 have {B1_bin_count}')
            print(f'bin belong to B2 have {B2_bin_count}')
            print(f'bin belong to B3 have {B3_bin_count}')
            print(f'bin belong to B4 have {B4_bin_count}')
            for compartment in compartment_ls:
                for i in range(-bin_radius, bin_radius+1):  # i is the relative distance to the center of the peak in bp
                    if compartment == 'A1':
                        csv_writer.writerow([experiment, compartment, i, A1_dic[i], A1_bin_count, line_number, 1000000 * A1_dic[i] / (A1_bin_count * line_number), ref_bins_name])
                    if compartment == 'A2':
                        csv_writer.writerow([experiment, compartment, i, A2_dic[i], A2_bin_count, line_number, 1000000 * A2_dic[i] / (A2_bin_count * line_number), ref_bins_name])
                    if compartment == 'B1':
                        csv_writer.writerow([experiment, compartment, i, B1_dic[i], B1_bin_count, line_number, 1000000 * B1_dic[i] / (B1_bin_count * line_number), ref_bins_name])
                    if compartment == 'B2':
                        csv_writer.writerow([experiment, compartment, i, B2_dic[i], B2_bin_count, line_number, 1000000 * B2_dic[i] / (B2_bin_count * line_number), ref_bins_name])
                    if compartment == 'B3':
                        csv_writer.writerow([experiment, compartment, i, B3_dic[i], B3_bin_count, line_number, 1000000 * B3_dic[i] / (B3_bin_count * line_number), ref_bins_name])
                    if compartment == 'B4':
                        csv_writer.writerow([experiment, compartment, i, B4_dic[i], B4_bin_count, line_number, 1000000 * B4_dic[i] / (B4_bin_count * line_number), ref_bins_name])
            # column names are:
            # experiment_name, compartment, position_to_nc_center, compartment_bin_count, compartment_bin_number
            # total_break_count, break_count/bin_count*total_break_count, reference_bin_name

    end_2 = time.perf_counter()
    print(f'{line_number} breaking position aggregated in {round(end_2 - start_2, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="to calculate aggreaged breaking frequecies relative to a buch of peak center positions")
    parser.add_argument("-b", help="peak signal centered indexed bin file", dest="ib", type=str, required=True)
    parser.add_argument("-r", help="bin radius", dest="radius", type=int, required=True)
    parser.add_argument("-p", help="breaking position file folder", dest="bpfd", type=str, required=True)
    parser.add_argument("-t", help="number of cpu cores", dest="t_number", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
