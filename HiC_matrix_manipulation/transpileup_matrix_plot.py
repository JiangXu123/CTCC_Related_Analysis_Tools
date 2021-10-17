#! /usr/bin/env python

import argparse
import numpy as np
import time
import os
import json
from matplotlib import pyplot as plt


def run(args):
    start_time = time.perf_counter()
    matrix_block_folder = args.mbf  # matrix_block_folder is the folder path that contain the blocked pileup matrices
    pad_size = args.pads
    random_trans_matrix_file = args.rtm
    TF_name_file = args.tf_n

    random_trans_matrix = np.load(random_trans_matrix_file)
    matrix_file_path_ls = []
    temp_matrix_file_name_ls = os.listdir(matrix_block_folder)
    for matrix_file in temp_matrix_file_name_ls:
        matrix_file_path_ls.append(os.path.abspath(matrix_block_folder) + '/' + matrix_file)

    A1_pileup_matrix = np.zeros((pad_size, pad_size), dtype=float)
    B1_pileup_matrix = np.zeros((pad_size, pad_size), dtype=float)
    B2_pileup_matrix = np.zeros((pad_size, pad_size), dtype=float)

    A1_matrix_count = 0
    B1_matrix_count = 0
    B2_matrix_count = 0

    average_A1_matrix = np.zeros((pad_size, pad_size), dtype=float)
    average_B1_matrix = np.zeros((pad_size, pad_size), dtype=float)
    average_B2_matrix = np.zeros((pad_size, pad_size), dtype=float)

    for matrix_path in matrix_file_path_ls:
        with open(matrix_path, 'r') as file1:
            matrix_obj = json.load(file1)
            A1_matrix_count += matrix_obj['A1'][1]
            B1_matrix_count += matrix_obj['B1'][1]
            B2_matrix_count += matrix_obj['B2'][1]
            A1_pileup_matrix += np.array(matrix_obj['A1'][0])
            B1_pileup_matrix += np.array(matrix_obj['B1'][0])
            B2_pileup_matrix += np.array(matrix_obj['B2'][0])

    if A1_matrix_count != 0:
        average_A1_matrix = A1_pileup_matrix / A1_matrix_count
    elif A1_matrix_count == 0:
        print('no valid coordinate found for A1-A1')
    if B1_matrix_count != 0:
        average_B1_matrix = B1_pileup_matrix / B1_matrix_count
    elif B1_matrix_count == 0:
        print('no valid coordinate found for B1-B1')
    if B2_matrix_count != 0:
        average_B2_matrix = B2_pileup_matrix / B2_matrix_count
    elif B2_matrix_count == 0:
        print('no valid coordinate found for B1-B1')

    label_coordinate_ls = [0, int((pad_size - 1) / 4), int((pad_size - 1) / 2), int(3 * (pad_size - 1) / 4), int((pad_size - 1))]
    label_ls = [-(pad_size - 1) / 2, -(pad_size - 1) / 4, 0, (pad_size - 1) / 4, (pad_size - 1) / 2]

    A1_subtracted_matrix = average_A1_matrix - random_trans_matrix
    B1_subtracted_matrix = average_B1_matrix - random_trans_matrix
    B2_subtracted_matrix = average_B2_matrix - random_trans_matrix

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(9, 3), dpi=300)
    color_map = plt.cm.get_cmap('PiYG')
    reversed_color_map = color_map.reversed()
    max_value = max([np.nanmax(A1_subtracted_matrix), np.nanmax(B1_subtracted_matrix), np.nanmax(B2_subtracted_matrix)])
    norm_thresholded = cm.colors.DivergingNorm(vmax=max_value, vcenter=0)

    axes[0].matshow(A1_subtracted_matrix, cmap=reversed_color_map, norm=norm_thresholded)
    axes[0].set_title(f"pileup over TF 1 on A1-A1 region", fontsize=8, y=1, pad=10)
    axes[0].set_xticks(label_coordinate_ls)
    axes[0].set_xticklabels(label_ls, size=6)
    axes[0].set_yticks(label_coordinate_ls)
    axes[0].set_yticklabels(label_ls, size=6)
    axes[0].tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False)

    axes[1].matshow(B1_subtracted_matrix, cmap=reversed_color_map, norm=norm_thresholded)
    axes[1].set_title(f"pileup over TF 1 on B1-B1 region", fontsize=8, y=1, pad=10)
    axes[1].set_xticks(label_coordinate_ls)
    axes[1].set_xticklabels(label_ls, size=6)
    axes[1].set_yticks(label_coordinate_ls)
    axes[1].set_yticklabels(label_ls, size=6)
    axes[1].tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False)

    axes[2].matshow(B2_subtracted_matrix, cmap=reversed_color_map, norm=norm_thresholded)
    axes[2].set_title(f"pileup over TF 1 on B2-B2 region", fontsize=8, y=1, pad=10)
    axes[2].set_xticks(label_coordinate_ls)
    axes[2].set_xticklabels(label_ls, size=6)
    axes[2].set_yticks(label_coordinate_ls)
    axes[2].set_yticklabels(label_ls, size=6)
    axes[2].tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False)
    end_time = time.perf_counter()

    print(f"matrix pileup combined and averaged in {round(end_time-start_time, 5)} seconds")


def main():
    parser = argparse.ArgumentParser(description="Pool and average the blocked pileup matrices, and subtract the random pileup matrix")
    parser.add_argument("-mf", help="blocked pileup matrix folder path", dest="mbf", type=str, required=True)
    parser.add_argument("-cs", help="chromosome size file", dest="chrsz", type=str, required=True)
    parser.add_argument("-rm", help="the random transpileup matrix", dest="rtm", type=str, required=True)
    parser.add_argument("-ps", help="pad size in terms of pixels", dest="pads", type=int, required=True)

    # parser.add_argument("-cn", help="the control background trans submatrix number", dest="rmn", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()