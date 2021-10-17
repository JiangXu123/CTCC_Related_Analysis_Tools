#! /usr/bin/env python

import argparse
import cooler
import tabix
import csv
import pandas as pd
import time
import subprocess
import concurrent.futures
import os


def run(args):
    start_time = time.perf_counter()
    cool_file = args.cool  # the cool file should be in the format .mcool::/resolutions
    TF_bed_file_1 = args.TF_bed_1  # the bed file contains 10 column, the 9th column contains the qValue (FDR), if the value is -1, then there's no FDR available
    TF_bed_file_2 = args.TF_bed_2  # if TF_bed_file_1 and TF_bed_file_2 are the same, only one TF mediated contact position will be considered
    compartment_file = args.comp_file  # this is an tabix indexed subcompartment file to categorize the ChIP signals
    TF_name_list_file = args.tf_l  # this list contains the TF number and the corresponding TF name, will be used to determine the TF name according to the number
    chrom_size_file = args.chrsz  # chromosome size file need to be editted to ensure only those chromosomes of interest are kept
    q_value_threshold = args.qth  # 1.3011 corresponding to FDR of 0.05
    pad_size = args.pad  # should be a odd number(the center pixel takes one pixel
    thread_number = args.t_number  # the thread_number will determine how many blocks the coordinate file be split to

    chromosome_ls = []
    with open(chrom_size_file, 'r') as file1:
        csv_reader_chr = csv.reader(file1, delimiter='\t')
        for line in csv_reader_chr:
            chromosome_ls.append(line[0])

    resolution = cool_file.split('/')[-1]
    tb = tabix.open(compartment_file)
    trans_comp_ls = ['A1', 'B1', 'B2']

    TF_no_name_dic = {}
    with open(TF_name_list_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        for line in csv_reader:
            TF_no_name_dic[line[0]] = line[1].split("_")[0]

    TF_1_name = TF_no_name_dic[TF_bed_file_1.split("/")[-1].split(".")[0]]
    TF_2_name = TF_no_name_dic[TF_bed_file_2.split("/")[-1].split(".")[0]]

    TF_1_bed_comp_file_name = TF_bed_file_1.split('.')[0] + '_' + 'comp' + '.bed'
    TF_2_bed_comp_file_name = TF_bed_file_2.split('.')[0] + '_' + 'comp' + '.bed'

    Trans_TF_looping_file_name = TF_bed_file_1.split('.')[0] + '_' + TF_bed_file_2.split('.')[0] + '.coord'
    Trans_TF_looping_folder_name = TF_bed_file_1.split('.')[0] + '_' + TF_bed_file_2.split('.')[0] + f'_{resolution}' + f'_{str(pad_size)}' # the folder_name looks like: 1_1_10000_201
    Trans_TF_looping_folder_path = os.path.abspath('./' + Trans_TF_looping_folder_name)
    Trans_TF_looping_file_path = os.path.abspath('./' + Trans_TF_looping_folder_name + '/' + Trans_TF_looping_file_name)
    try:
        os.makedirs(Trans_TF_looping_folder_path)  # create a trans_looping folder to contain the coordinate and result
        print(f"Folder for TF looping calculation created at {Trans_TF_looping_folder_name}")
    except FileExistsError:
        print(f'Folder for TF looping calculation already exists at {Trans_TF_looping_folder_name}, will use for calculation')
    splitting_coordinate_temp_folder_path = Trans_TF_looping_folder_path + '/split'
    os.mkdir(splitting_coordinate_temp_folder_path)

    TF_1_comp_bed_file_path = os.path.abspath('./' + Trans_TF_looping_folder_name + '/' + TF_1_bed_comp_file_name)
    TF_2_comp_bed_file_path = os.path.abspath('./' + Trans_TF_looping_folder_name + '/' + TF_2_bed_comp_file_name)

    with open(TF_bed_file_1, 'r') as file1:
        with open(TF_1_comp_bed_file_path, 'w') as file2:
            csv_reader_1 = csv.reader(file1, delimiter='\t')
            csv_writer_1 = csv.writer(file2, delimiter='\t')
            first_line_qvalue_1 = float(next(csv_reader_1)[8])
            if first_line_qvalue_1 == -1:  # if q value is not valid, all position will be taken into account
                for line in csv_reader_1:
                    try:
                        query_results = tb.query(line[0], int(line[1]), int(line[2]))
                        result = next(query_results)
                        for comp_1 in trans_comp_ls:
                            if result[3] == comp_1:
                                csv_writer_1.writerow([line[0], int(line[1]), int(line[2]), comp_1])
                    except:
                        pass
            elif first_line_qvalue_1 != -1:  # if q value is valid, then position with qvalue > q_value_threshold will be considered
                for line in csv_reader_1:
                    if float(line[8]) > q_value_threshold:  # coresponding q value < 0.05
                        try:
                            query_results = tb.query(line[0], int(line[1]), int(line[2]))
                            result = next(query_results)
                            for comp_2 in trans_comp_ls:
                                if result[3] == comp_2:
                                    csv_writer_1.writerow([line[0], int(line[1]), int(line[2]), comp_2])
                        except:
                            pass

    with open(TF_bed_file_2, 'r') as file3:
        with open(TF_2_comp_bed_file_path, 'w') as file4:
            csv_reader_2 = csv.reader(file3, delimiter='\t')
            csv_writer_2 = csv.writer(file4, delimiter='\t')
            first_line_qvalue_2 = float(next(csv_reader_2)[8])
            if first_line_qvalue_2 == -1:  # if q value is not valid, all position will be taken into account
                for line in csv_reader_2:
                    try:
                        query_results = tb.query(line[0], int(line[1]), int(line[2]))
                        result = next(query_results)
                        for comp_3 in trans_comp_ls:
                            if result[3] == comp_3:
                                csv_writer_2.writerow([line[0], int(line[1]), int(line[2]), comp_3])
                    except:
                        pass
            elif first_line_qvalue_2 != -1:  # if q value is valid, then position with qvalue > q_value_threshold will be considered
                for line in csv_reader_2:
                    if float(line[8]) > q_value_threshold:  # coresponding q value < 0.05
                        try:
                            query_results = tb.query(line[0], int(line[1]), int(line[2]))
                            result = next(query_results)
                            for comp_4 in trans_comp_ls:
                                if result[3] == comp_4:
                                    csv_writer_2.writerow([line[0], int(line[1]), int(line[2]), comp_4])
                        except:
                            pass

    TF_1_df = pd.read_csv(TF_1_comp_bed_file_path, delimiter='\t', names=['chr_name', 'start', 'end', 'subcompartment'])
    TF_2_df = pd.read_csv(TF_2_comp_bed_file_path, delimiter='\t', names=['chr_name', 'start', 'end', 'subcompartment'])
    '''generate a trans interaction list that is based on the same-subcompartment(A1-A1, B1-B1, B2-B2) and two chip-seq signal file, the positions in this list
      do not consider the size fo the submatrix, and will be filtered in later step'''
    with open(Trans_TF_looping_file_path, 'w') as file1:
        csv_writer = csv.writer(file1, delimiter='\t')
        for comp in trans_comp_ls:
            comp_TF_1_df = TF_1_df.loc[TF_1_df['subcompartment'] == comp]
            comp_TF_2_df = TF_2_df.loc[TF_2_df['subcompartment'] == comp]
            for i in range(0, len(chromosome_ls)):
                for j in range(i + 1, len(chromosome_ls)):
                    intermediate_time_1 = time.perf_counter()
                    chr_1_name = chromosome_ls[i]
                    chr_1_positions_df_1 = comp_TF_1_df.loc[comp_TF_1_df['chr_name'] == chr_1_name]
                    chr_2_name = chromosome_ls[j]
                    chr_2_positions_df_2 = comp_TF_2_df.loc[comp_TF_2_df['chr_name'] == chr_2_name]
                    for index_1, row_1 in chr_1_positions_df_1.iterrows():
                        for index_2, row_2 in chr_2_positions_df_2.iterrows():
                            csv_writer.writerow([row_1[0], row_1[1], row_1[2], row_2[0], row_2[1], row_2[2], comp])
                    intermediate_time_2 = time.perf_counter()
                    print(f"interaction list for chromosome {chromosome_ls[i]} against {chromosome_ls[j]} between subcompartment {comp} in {round(intermediate_time_2 - intermediate_time_1, 5)} seconds")
    end_time = time.perf_counter()
    print(f"trans interaction list for {TF_1_name} and {TF_2_name} with FDR value less than {10 ** (-q_value_threshold)}have been created in {round(end_time - start_time, 5)} seconds")

    get_total_coordinate_num = subprocess.check_output(f'cat {Trans_TF_looping_file_path} | wc -l', shell=True)
    total_coordinate_num = int(get_total_coordinate_num.split()[0].decode("utf-8"))  # total_coordinate_num is the total number of coordinates representing the trans same compartment contact
    block_lines = int(total_coordinate_num / thread_number)  # make an estimation how many lines should contain per block
    print(f'total of {total_coordinate_num} coordinates will be split into {block_lines} per temporary block file')

    os.chdir(splitting_coordinate_temp_folder_path)  # enter the temporary folder for splitting the coordinate file
    subprocess.call(f'cat {Trans_TF_looping_file_path} | split -d -l {block_lines} - temp_in_', shell=True)  # I cheated here using bash shell's split method, which is much faster in splitting files
    split_block_num = len(os.listdir(splitting_coordinate_temp_folder_path))
    print(f'coordinate file  {Trans_TF_looping_file_name} is split to {split_block_num} blocks into the temp folder')

    temp_block_file_name_ls = os.listdir(splitting_coordinate_temp_folder_path)
    temp_block_file_path_ls = []

    for name in temp_block_file_name_ls:
        temp_block_file_path_ls.append(splitting_coordinate_temp_folder_path + '/' + name)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        processes = [executor.submit(trans_coordinate_generator, cool_file, coordinate_block_file_path, pad_size, chromosome_ls) for coordinate_block_file_path in temp_block_file_path_ls]
        for process in concurrent.futures.as_completed(processes):
            print(f'for {process.result()[0]}, {process.result()[1]} trans coordinates were parsed, finding {process.result()[2]} valid coordinate that fit the trans pileup condition，taking {process.result()[3]} seconds')
    for input_block_file in temp_block_file_path_ls:
        subprocess.call(f'rm -rf {input_block_file}', shell=True)

    # valid_trans_coordinate_block_ls = os.listdir(splitting_coordinate_temp_folder_path)
    # valid_trans_coordinate_block_path_ls = []
    # for valid_coord_block_name in valid_trans_coordinate_block_ls:
    #     valid_trans_coordinate_block_path_ls.append(splitting_coordinate_temp_folder_path + '/' + valid_coord_block_name)

    # with open(Valid_trans_coord_file_path, 'a') as file1:  # concatenate the valid_trans_coordinate_block_file to a single file
    #     csv_writer = csv.writer(file1, delimiter='\t')
    #     for valid_trans_coordinate_block_path in valid_trans_coordinate_block_path_ls:
    #         with open(valid_trans_coordinate_block_path, 'r') as file2:
    #             csv_reader = csv.reader(file2, delimiter='\t')
    #             for line in csv_reader:
    #                 csv_writer.writerow(line)
    # shutil.rmtree(splitting_coordinate_temp_folder_path)  # remove the temp folder for coordinate splitting and processing, including the blocked input and output file
    # subprocess.call(f'rm -rf {Trans_TF_looping_file_path}', shell=True) # remove the unfiltered Trans_TF_looping_file
    end_time = time.perf_counter()
    print(f'process finished, takes {round(end_time-start_time, 5)} seconds')


def trans_coordinate_generator(cool_file, coordinate_file, pad, chromosome_ls):
    process_time_start = time.perf_counter()
    coordinate_file_name = coordinate_file.split('/')[-1]
    qualified_coordinate_file_file_path = coordinate_file + '.' + 'qualified'
    c1 = cooler.Cooler(cool_file)
    bins = c1.bins()[:]
    chr_bin_dic = {}

    for chromosome in chromosome_ls:
        chr_bin_dic[chromosome] = (bins.loc[bins['chrom'] == chromosome].index[0], bins.loc[bins['chrom'] == chromosome].index[-1])  # generate a tuple to store the chromosome's start and end bin number
    line_count = 0
    valid_coord_count = 0
    with open(coordinate_file, 'r') as file1:
        with open(qualified_coordinate_file_file_path, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            for line in csv_reader:
                mat_center_x_bin_no = -1
                mat_center_y_bin_no = -1
                line_count += 1
                chr_1_name = line[0]
                chr_1_start = int(line[1])
                chr_1_end = int(line[2])

                mat_cent_x_pos = int(chr_1_start + (chr_1_end - chr_1_start) / 2)
                try:
                    mat_center_x_bin = bins.loc[(bins['chrom'] == chr_1_name) & (mat_cent_x_pos > bins['start']) & (mat_cent_x_pos < bins['end'])]
                    mat_center_x_bin_no = mat_center_x_bin.index[0]
                except IndexError:
                    pass
                chr_2_name = line[3]
                chr_2_start = int(line[4])
                chr_2_end = int(line[5])

                mat_cent_y_pos = int(chr_2_start + (chr_2_end - chr_2_start) / 2)
                try:
                    mat_center_y_bin = bins.loc[(bins['chrom'] == chr_2_name) & (mat_cent_y_pos > bins['start']) & (mat_cent_y_pos < bins['end'])]
                    mat_center_y_bin_no = mat_center_y_bin.index[0]
                except IndexError:
                    pass
                if (mat_center_x_bin_no != -1) & (mat_center_y_bin_no != -1):  # if a coordinate is found in the both the x and y bin table
                    matrix_r = int((pad - 1) / 2)  # matrix_r is the pad half width not including the central pixel.
                    mat_x_start_bin_no = int(mat_center_x_bin_no - matrix_r)
                    mat_x_end_bin_no = int(mat_center_x_bin_no + matrix_r)
                    mat_y_start_bin_no = int(mat_center_y_bin_no - matrix_r)
                    mat_y_end_bin_no = int(mat_center_y_bin_no + matrix_r)
                    if (mat_x_start_bin_no >= chr_bin_dic[chr_1_name][0]) & (mat_x_end_bin_no <= chr_bin_dic[chr_1_name][1]) & (mat_y_start_bin_no >= chr_bin_dic[chr_2_name][0]) & (mat_y_end_bin_no <= chr_bin_dic[chr_2_name][1]):
                        matrix_x_chr = chr_1_name
                        matrix_x_start = bins.iloc[mat_x_start_bin_no]['start']
                        matrix_x_end = bins.iloc[mat_x_end_bin_no]['end']
                        matrix_y_chr = chr_2_name
                        matrix_y_start = bins.iloc[mat_y_start_bin_no]['start']
                        matrix_y_end = bins.iloc[mat_y_end_bin_no]['end']
                        csv_writer.writerow([matrix_x_chr, matrix_x_start, matrix_x_end, matrix_y_chr, matrix_y_start, matrix_y_end, line[6]])
                        valid_coord_count += 1
            process_time_end = time.perf_counter()
            total_process_time = process_time_end - process_time_start
            return coordinate_file_name, line_count, valid_coord_count, total_process_time


def main():
    parser = argparse.ArgumentParser(description="To do pileups with ChIP narrow bed file in the trans region")
    parser.add_argument("-c", help="cool file to be analyzed", dest="cool", type=str, required=True)
    parser.add_argument("-b1", help="narrow bed file of the chip-seq data for TF1", dest="TF_bed_1", type=str, required=True)
    parser.add_argument("-b2", help="narrow bed file of the chip-seq data for TF2", dest="TF_bed_2", type=str, required=True)
    parser.add_argument("-cp", help="subcompartment file", dest="comp_file", type=str, required=True)
    parser.add_argument("-nl", help="TF number name list file", dest="tf_l", type=str, required=True)
    parser.add_argument("-cs", help="chromosome size file", dest="chrsz", type=str, required=True)
    parser.add_argument("-qt", help="q value threshold for qualified ChIP peaks, suggestion: 1.3011 corresponding to FDR of 0.05", dest="qth", type=float, required=True)
    parser.add_argument("-pd", help="pad size(submatrix diemention (square matrix: pad*pad) in terms of pixels", dest="pad", type=int, required=True)
    parser.add_argument("-tn", help="number of cpu core used, process(thread) number", dest="t_number", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
