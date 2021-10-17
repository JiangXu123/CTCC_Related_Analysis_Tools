#! /usr/bin/env python

import argparse
import cooler
import numpy as np
import csv
import time
import subprocess
import shutil
import concurrent.futures
import os
import json


def run(args):
    start_time = time.perf_counter()
    valid_trans_coordinate_file = args.vcoord
    cool_file = args.cool  # the cool file should be in the format .mcool::/resolutions
    thread_number = args.t_number  # the thread_number will determine how many blocks the coordinate file be split to
    chrom_size_file = args.chrsz
    pad_size = args.pads
    chromosome_ls = []
    with open(chrom_size_file, 'r') as file1:
        csv_reader_chr = csv.reader(file1, delimiter='\t')
        for line in csv_reader_chr:
            chromosome_ls.append(line[0])
    pileup_info_name = valid_trans_coordinate_file.split('.')[0]
    temporary_folder_name = pileup_info_name + '_temp'
    os.mkdir(temporary_folder_name)  # each temporary_folder has their unique name and will be deleted once pileup work finish
    valid_coordinate_folder_path = os.path.abspath('.')
    temporary_folder_path = valid_coordinate_folder_path + '/' + temporary_folder_name
    valid_trans_coordinate_file_path =  valid_coordinate_folder_path + '/' + valid_trans_coordinate_file
    pileup_matrix_folder = valid_coordinate_folder_path + '/pileup_matrices'
    try:
        os.mkdir(pileup_matrix_folder)  # the pileup_matrix_folder is shared with different TF, pad and resolutions
    except FileExistsError:
        print('pileup_matrices folder already exists, will use to store new pileup matrices')
        pass
    get_total_coordinate_num = subprocess.check_output(f'cat {valid_trans_coordinate_file_path} | wc -l', shell=True)
    total_valid_coordinate_num = int(get_total_coordinate_num.split()[0].decode("utf-8"))  # total_coordinate_num is the total number of coordinates representing the trans same compartment contact
    block_lines = int(total_valid_coordinate_num / thread_number)  #  make an estimation how many lines should contain per block
    os.chdir(temporary_folder_path)  # enter the temporary folder for splitting the valid_coordinate file
    subprocess.call(f'cat {valid_trans_coordinate_file_path} | split -d -l {block_lines} - temp_valid_coord_block_', shell=True)
    split_block_num = len(os.listdir(temporary_folder_path))
    print(f'valid_coordinate file  {valid_trans_coordinate_file} is split to {split_block_num} blocks into the temp folder {temporary_folder_name}')

    temp_block_file_name_ls = os.listdir(temporary_folder_path)
    temp_block_file_path_ls = []
    for name in temp_block_file_name_ls:
        temp_block_file_path_ls.append(temporary_folder_path + '/' + name)

    pile_up_matrix = {'A1': np.zeros((pad_size, pad_size)), 'B1': np.zeros((pad_size, pad_size)), 'B2': np.zeros((pad_size, pad_size))}
    pile_up_matrix_count = {'A1': 0, 'B1': 0, 'B2': 0}

    with concurrent.futures.ProcessPoolExecutor() as executor:
        processes = [executor.submit(matrix_pileup, cool_file, coordinate_block_file_path, pad_size, pileup_info_name) for coordinate_block_file_path in temp_block_file_path_ls]
        for process in concurrent.futures.as_completed(processes):
            pile_up_matrix['A1'] += process.result()[0]['A1']  # pileup the result from each matrix to the
            pile_up_matrix['B1'] += process.result()[0]['B1']
            pile_up_matrix['B2'] += process.result()[0]['B2']
            pile_up_matrix_count['A1'] += process.result()[1]['A1']
            pile_up_matrix_count['B1'] += process.result()[1]['B1']
            pile_up_matrix_count['B2'] += process.result()[1]['B2']
            if np.isnan(process.result()[0]['A1']).any():  # 0 divide 0 will result in numpy.nan, any integer divide 0 will result in numpy.inf
                print('nan found in pileup matrix A1')
            if np.isnan(process.result()[0]['B1']).any():  # a.any() means if any numpy.nan in the numpy array
                print('nan found in pileup matrix B1')
            if np.isnan(process.result()[0]['B2']).any():
                print('nan found in pileup matrix B2')
    sub_pile_up_matrix = {'A1': (pile_up_matrix['A1'], pile_up_matrix_count['A1']), 'B1': (pile_up_matrix['B1'], pile_up_matrix_count['B1']), 'B2': (pile_up_matrix['B2'], pile_up_matrix_count['B2'])}
    shutil.rmtree(temporary_folder_path) # remove the temporary folder for splitting the matrix
    os.chdir(pileup_matrix_folder)
    saved_pileup_dic_name = pileup_info_name + '.json'

    class NumpyEncoder(json.JSONEncoder):   # to enable json handle numpy array
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            return json.JSONEncoder.default(self, obj)

    with open(saved_pileup_dic_name, 'w') as file_obj:
        json.dump(sub_pile_up_matrix, file_obj, cls=NumpyEncoder) # pileup will be stored in json format, containing subcompartment information, pileup matrices sum, and pileup matrices count

    end_time = time.perf_counter()
    print(f'Trans matrix pileup generate {total_valid_coordinate_num} submatrices, with {thread_number} cores of CPU, taking {round(end_time-start_time, 5)} seconds')


def matrix_pileup(cool_file, valid_coordinate_file, pad, pileup_name):
    sub_matrix_pileup = {'A1': np.zeros((pad, pad), dtype='float'), 'B1': np.zeros((pad, pad), dtype='float'), 'B2': np.zeros((pad, pad), dtype='float')}
    sub_matrix_count = {'A1':0, 'B1':0, 'B2':0}
    c1 = cooler.Cooler(cool_file)
    cool_file_name = cool_file.split('.')[0]
    get_total_coordinate_num = subprocess.check_output(f'cat {valid_coordinate_file} | wc -l', shell=True)
    block_coordinate_num = int(get_total_coordinate_num.split()[0].decode("utf-8"))
    pileup_time_1 = time.perf_counter()
    resolution = cool_file.split('/')[-1]
    line_count = 0
    print(f'start pileup for {valid_coordinate_file}')
    with open(valid_coordinate_file, 'r') as file1:
        csv_reader = csv.reader(file1, delimiter='\t')
        for line in csv_reader:
            chr_x_name = line[0]
            chr_x_start = int(line[1])
            chr_x_end = int(line[2])
            chr_y_name = line[3]
            chr_y_start = int(line[4])
            chr_y_end = int(line[5])
            subcompartment_type = line[6]
            sub_matrix = c1.matrix(balance=True, sparse=False).fetch((chr_x_name, chr_x_start, chr_x_end), (chr_y_name, chr_y_start, chr_y_end))
            sub_matrix_pileup[subcompartment_type] += np.nan_to_num(sub_matrix)
            sub_matrix_count[subcompartment_type] += 1
            line_count += 1
            if line_count%1000 == 0:
                pileup_time_2 = time.perf_counter()
                print(f'for {pileup_name} with {cool_file_name} with {pad} pad, {resolution} pixel length, from split file {valid_coordinate_file} with {block_coordinate_num} valid coordinates, {line_count} has been parsed, taking {round(pileup_time_2-pileup_time_1, 5)} seconds passed')
    return sub_matrix_pileup, sub_matrix_count


def main():
    parser = argparse.ArgumentParser(description="To do pileups in the trans region")
    parser.add_argument("-c", help="cool file to be analyzed", dest="cool", type=str, required=True)
    parser.add_argument("-vc", help="valid coordinate file for trans pileup", dest="vcoord", type=str, required=True)
    parser.add_argument("-cs", help="chromosome size file", dest="chrsz", type=str, required=True)
    parser.add_argument("-tn", help="number of cpu core used, process(thread) number", dest="t_number", type=int, required=True)
    parser.add_argument("-ps", help="pad size in terms of pixels", dest="pads", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()