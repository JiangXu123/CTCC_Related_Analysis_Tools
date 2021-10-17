#! /usr/bin/env python

import csv
import argparse
import cooler
import numpy as np
import time
import random
import os


def run(args):
    start_time = time.perf_counter()
    cool_file = args.cool  # the cool file should be in the format .mcool::/resolutions, and should contain the absolute path
    cool_file_name = cool_file.split('.')[0].split('/')[-1]  # the cool file
    rand_mat_num = args.rmn  # This number is assigned arbitrarily, will be improved in future based on using experience
    chrom_size_file = args.chrsz
    pad_size = args.pads
    resolution = cool_file.split('/')[-1]
    random_trans_matrix_pileup_name = cool_file_name + '_' + str(resolution) + '_' + str(pad_size) + '_' + str(rand_mat_num) + '_' + 'trans_background'
    # a trans background matrix name should look like: 4DNFIXP4QG5B_10000_201_3000_trans_background.npy
    # here the first number is the resolution, the second is pad and the third is the random number of matrices used
    chromosome_ls = []
    with open(chrom_size_file, 'r') as file1:
        csv_reader_chr = csv.reader(file1, delimiter='\t')
        for line in csv_reader_chr:
            chromosome_ls.append(line[0])
    trans_background_pileup_matrix = random_trans_matrix_generator(rand_mat_num, pad_size, cool_file, chromosome_ls)
    trans_background_pileup_matrix_file_path = os.path.abspath('.') + '/' + random_trans_matrix_pileup_name
    np.save(trans_background_pileup_matrix_file_path, trans_background_pileup_matrix)


def random_trans_matrix_generator(random_matrix_num, pad, cool_file, chromosome_ls): # random_matrix_num is the number of random trans matrix needed for background calculation
    c1 = cooler.Cooler(cool_file)
    # chrom_ls_from_cool_file = c1.chromsizes.reset_index()['name'].to_list
    # filtered_out_chrom_ls = list(set(chrom_ls_from_cool_file)-set(chromosome_ls))
    bins = c1.bins()[:]
    GM12878_bins = bins.loc[(bins['chrom'] != 'chrY') & (bins['chrom'] != 'chrM')]
    chr_bin_dic = {}
    for chromosome in chromosome_ls:
        chr_bin_dic[chromosome] = (bins.loc[bins['chrom'] == chromosome].index[0], bins.loc[bins['chrom'] == chromosome].index[-1]) # chr_bin_dic looks like {chr1:(1,280878789), chr2:...}
    total_bin_num = len(GM12878_bins.index)
    ctrl_mat_bin_ls = []
    while len(ctrl_mat_bin_ls) < random_matrix_num:
        potential_rando_x_bin = random.randint(0, total_bin_num)
        potential_rando_y_bin = random.randint(0, total_bin_num)
        poteintial_rando_x_bin_chr = GM12878_bins.iloc[potential_rando_x_bin]['chrom']
        poteintial_rando_y_bin_chr = GM12878_bins.iloc[potential_rando_y_bin]['chrom']
        poteintial_rando_x_bin_chr_bin_start = chr_bin_dic[poteintial_rando_x_bin_chr][0]
        poteintial_rando_x_bin_chr_bin_end = chr_bin_dic[poteintial_rando_x_bin_chr][1]
        poteintial_rando_y_bin_chr_bin_start = chr_bin_dic[poteintial_rando_y_bin_chr][0]
        poteintial_rando_y_bin_chr_bin_end = chr_bin_dic[poteintial_rando_y_bin_chr][1]
        if (poteintial_rando_x_bin_chr != poteintial_rando_y_bin_chr) & ((potential_rando_x_bin - ((pad - 1) / 2)) >= poteintial_rando_x_bin_chr_bin_start) & ((potential_rando_x_bin + ((pad - 1) / 2)) <= poteintial_rando_x_bin_chr_bin_end) & ((potential_rando_y_bin - ((pad - 1) / 2)) >= poteintial_rando_y_bin_chr_bin_start) & ((potential_rando_y_bin + ((pad - 1) / 2)) <= poteintial_rando_y_bin_chr_bin_end):
            ctrl_mat_bin_ls.append((potential_rando_x_bin, potential_rando_y_bin))
    print(f"{len(ctrl_mat_bin_ls)} trans random matrix coordinate generated as set number to {random_matrix_num}")
    ctrl_pile_up_matrix = np.zeros((pad, pad), dtype=float)
    actual_rand_trans_matrix_count = 0
    for coordinate in ctrl_mat_bin_ls:  # coordinate is a tuple (x, y) in terms of bin's number
        ctrl_x_chr = bins.iloc[coordinate[0]]['chrom']
        ctrl_x_chr_start = bins.iloc[int(coordinate[0] - ((pad-1)/2))]['start']
        ctrl_x_chr_end = bins.iloc[int(coordinate[0] + ((pad-1)/2))]['end']
        ctrl_y_chr = bins.iloc[coordinate[1]]['chrom']
        ctrl_y_chr_start = bins.iloc[int(coordinate[1] - ((pad - 1)/2))]['start']
        ctrl_y_chr_end = bins.iloc[int(coordinate[1] + ((pad - 1)/2))]['end']
        ctrl_sub_matrix = c1.matrix(balance=True, sparse=False).fetch((ctrl_x_chr, ctrl_x_chr_start, ctrl_x_chr_end), (ctrl_y_chr, ctrl_y_chr_start, ctrl_y_chr_end))
        ctrl_pile_up_matrix += np.nan_to_num(ctrl_sub_matrix)
        actual_rand_trans_matrix_count += 1
    return ctrl_pile_up_matrix/actual_rand_trans_matrix_count


def main():
    parser = argparse.ArgumentParser(description="To generate background pileups at trans pileup at random locations")
    parser.add_argument("-c", help="cool file to be analyzed", dest="cool", type=str, required=True)
    parser.add_argument("-cs", help="chromosome size file", dest="chrsz", type=str, required=True)
    parser.add_argument("-ps", help="pad size in terms of pixels", dest="pads", type=int, required=True)
    parser.add_argument("-cn", help="the control background trans submatrix number", dest="rmn", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()