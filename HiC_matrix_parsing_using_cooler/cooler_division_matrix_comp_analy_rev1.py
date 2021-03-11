#! /usr/bin/env python

import argparse
import cooler
import tabix
import numpy as np
import csv
import pandas as pd
import gzip
import time


def run(args):
    cool_file_1 = args.cool_1
    cool_file_2 = args.cool_2
    division_cis_pixel_value_file = args.div_cis_pix
    division_trans_pixel_value_file = args.div_trans_pix
    div_matrix_pixel_file = args.div_matrix_thresholded_pix
    compartment_file = args.comp_file
    increment_step = args.incre_step
    aggr_comp_pixel_count = args.aggr_th_comp_ct
    pixel_length_1 = int(cool_file_1.split("/")[-1])
    pixel_length_2 = int(cool_file_2.split("/")[-1])
    pixel_length = pixel_length_1
    if pixel_length_1 != pixel_length_2:
        print(f"input cool file resolution is different, please check input file!")
    tb = tabix.open(compartment_file)
    c1 = cooler.Cooler(cool_file_1)
    c2 = cooler.Cooler(cool_file_2)
    start = time.perf_counter()
    # generate a chromosome list to contain all the chromosomes
    chromosome_ls = []
    for i in range(1, 23):
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')
    chromosome_ls.append('chrY')
    chromosome_ls.append('chrM')
    # the chromosome list is to contain the chromosomes that will be investigated
    # first step is to calculate the cis and trans average value(devoid of nan, inf and 0) of the division matrix

    with open(division_cis_pixel_value_file, "w") as f1:
        csv_writer = csv.writer(f1, delimiter='\t')
        csv_writer.writerow(["chr_name", "value"])
        for i in range(0, len(chromosome_ls)):
            cis_1 = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[i])
            cis_2 = c2.matrix(sparse=False, balance=True).fetch(chromosome_ls[i])
            cis_div_matrix = np.divide(cis_1, cis_2)
            for value in np.nditer(cis_div_matrix):
                if (not np.isinf(value)) & (not np.isnan(value)) & (value != 0):  # get rid of 0, inf, and nan to avoid the disturbance from these values
                    csv_writer.writerow([chromosome_ls[i], value])

    with open(division_trans_pixel_value_file, "w") as f2:
        csv_writer = csv.writer(f2, delimiter='\t')
        csv_writer.writerow(["chr_1_name", "chr_2_name", "value"])
        for i in range(0, len(chromosome_ls)):
            for j in range(i + 1, len(chromosome_ls)):
                trans_1 = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[1], chromosome_ls[2])
                trans_2 = c2.matrix(sparse=False, balance=True).fetch(chromosome_ls[1], chromosome_ls[2])
                trans_div_matrix = np.divide(trans_1, trans_2)
                for value in np.nditer(trans_div_matrix):
                    if (not np.isinf(value)) & (not np.isnan(value)) & (value != 0):  # get rid of 0, inf, and nan to avoid the disturbance from these values
                        csv_writer.writerow([chromosome_ls[i], chromosome_ls[j], value])

    cis_pixel_df = pd.read_csv(division_cis_pixel_value_file, delimiter='\t')
    max_cis_pixel_value = cis_pixel_df.describe().loc['max']['value']
    min_cis_pixel_value = cis_pixel_df.describe().loc['min']['value']

    trans_pixel_df = pd.read_csv(division_trans_pixel_value_file, delimiter='\t')
    max_trans_pixel_value = trans_pixel_df.describe().loc['max']['value']
    min_trans_pixel_value = trans_pixel_df.describe().loc['min']['value']
    print(f"The highest cis value of the division matrix is {max_cis_pixel_value}")
    print(f"The lowest cis value of the division matrix is {min_cis_pixel_value}")
    print(f"The highest trans value of the division matrix is {max_trans_pixel_value}")
    print(f"The lowest trans value of the division matrix is {min_trans_pixel_value}")

    # comp_len_dic contains the length of subcompartments (bp) of the GM12878 cell line
    # cis_pixel_total_counts_dic contains the cis same compartment pixel counts for each subcompartment,
    # used to normalize the data considering the difference in compartments' length
    comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}

    with gzip.open(compartment_file, 'rt') as f3:
        csv_reader = csv.reader(f3, delimiter='\t')
        chr_comp_dic = {}
        for chromosome in chromosome_ls:
            chr_comp_dic[chromosome] = {}
            chr_comp_dic[chromosome]['A1'] = 0
            chr_comp_dic[chromosome]['A2'] = 0
            chr_comp_dic[chromosome]['B1'] = 0
            chr_comp_dic[chromosome]['B2'] = 0
            chr_comp_dic[chromosome]['B3'] = 0
            chr_comp_dic[chromosome]['B4'] = 0
            chr_comp_dic[chromosome]['NA'] = 0
        for line in csv_reader:
            chr_comp_dic[line[0]][line[3]] += (int(line[2]) - int(line[1]))  # generate each chomosome's compartment length
            comp_len_dic[line[3]] += (int(line[2]) - int(line[1]))  # generate each compartment length

    cis_pixel_total_counts_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0, 'NA': 0}
    for chromosome in chromosome_ls:
        cis_pixel_total_counts_dic['A1'] += int(chr_comp_dic[chromosome]['A1'] / pixel_length) * int(chr_comp_dic[chromosome]['A1'] / pixel_length)
        cis_pixel_total_counts_dic['A2'] += int(chr_comp_dic[chromosome]['A2'] / pixel_length) * int(chr_comp_dic[chromosome]['A2'] / pixel_length)
        cis_pixel_total_counts_dic['B1'] += int(chr_comp_dic[chromosome]['B1'] / pixel_length) * int(chr_comp_dic[chromosome]['B1'] / pixel_length)
        cis_pixel_total_counts_dic['B2'] += int(chr_comp_dic[chromosome]['B2'] / pixel_length) * int(chr_comp_dic[chromosome]['B2'] / pixel_length)
        cis_pixel_total_counts_dic['B3'] += int(chr_comp_dic[chromosome]['B3'] / pixel_length) * int(chr_comp_dic[chromosome]['B3'] / pixel_length)
        cis_pixel_total_counts_dic['B4'] += int(chr_comp_dic[chromosome]['B4'] / pixel_length) * int(chr_comp_dic[chromosome]['B4'] / pixel_length)
        cis_pixel_total_counts_dic['NA'] += int(chr_comp_dic[chromosome]['NA'] / pixel_length) * int(chr_comp_dic[chromosome]['NA'] / pixel_length)

    total_comp_pixel_number_dic = {}  # total_comp_pixel_number_dic = trans_pixel_total_counts_dic + cis_pixel_total_counts_dic
    trans_pixel_total_counts_dic = {}
    comp_ls = ['A1', 'A2', 'B1', 'B2', 'B3', 'B4', 'NA']

    # to calculate the total pixel number belong to the same compartment interaction matrix, which include trans and cis
    for compartment in comp_ls:
        total_comp_pixel_number_dic[compartment] = int(comp_len_dic[compartment] / pixel_length) * int(comp_len_dic[compartment] / pixel_length)
        trans_pixel_total_counts_dic[compartment] = total_comp_pixel_number_dic[compartment] - cis_pixel_total_counts_dic[compartment]

    # to get the bin and corresponding locations as a dictionary
    bins = c1.bins()[:]  # get the bin dataframe from one of the Hic matrix
    bin_dic = {}
    for index, row in bins.iterrows():  # bins is a pandas dataframe containing the bin information from the cooler file
        bin_dic[index] = [row[0], row[1], row[2]]

    with open(div_matrix_pixel_file, "w") as f4:
        csv_writer = csv.writer(f4, delimiter='\t')
        csv_writer.writerow(['chr_1_name', 'chr_2_name', 'subcompartment', 'threshold', 'interaction_type'])  # write the header of the output file
        for k in range(0, len(chromosome_ls)):
            inter_time_start = time.perf_counter()
            cis_1 = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[k])
            cis_2 = c2.matrix(sparse=False, balance=True).fetch(chromosome_ls[k])
            cis_div_matrix = np.divide(cis_1,cis_2)
            for y in range(0, cis_div_matrix.shape[0]):
                for x in range(0, cis_div_matrix.shape[1]):
                    cis_pixel_value = cis_div_matrix[y, x]  # y,x are the row and column of the pixel in the individual trans matrix
                    try:
                        y_comp = next(tb.query(bin_dic[y][0], bin_dic[y][1], bin_dic[y][2]))[3]
                    except:
                        y_comp = "y_comp_nv"
                    try:
                        x_comp = next(tb.query(bin_dic[x][0], bin_dic[x][1], bin_dic[x][2]))[3]
                    except:
                        x_comp = "x_comp_nv"
                    if y_comp == x_comp:
                        for cis_pixel_value_threshold in np.arange(min_cis_pixel_value, max_cis_pixel_value, increment_step):
                            if cis_pixel_value >= cis_pixel_value_threshold:
                                csv_writer.writerow([chromosome_ls[k], chromosome_ls[k], y_comp, cis_pixel_value_threshold, 'cis'])
            inter_time_end = time.perf_counter()
            print(f"{chromosome_ls[k]} cis matrix parsed in {round(inter_time_end-inter_time_start, 2)} seconds")
        for i in range(0, len(chromosome_ls)): # generate the over-the-threshold trans pixel contents first
            for j in range(i + 1, len(chromosome_ls)):
                inter_time_start = time.perf_counter()
                trans_1 = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[i], chromosome_ls[j])
                trans_2 = c2.matrix(sparse=False, balance=True).fetch(chromosome_ls[i], chromosome_ls[j])
                trans_div_matrix = np.divide(trans_1, trans_2)
                chr_name_1 = chromosome_ls[i]  # first chromsome name
                chr_name_2 = chromosome_ls[j]  # second chromosome name
                chr_bins_1 = bins.loc[bins['chrom'] == chr_name_1]
                chr_bins_2 = bins.loc[bins['chrom'] == chr_name_2]
                chr_1_bin_abs_start = chr_bins_1.index[0]
                chr_2_bin_abs_start = chr_bins_2.index[0]
                for y in range(0, trans_div_matrix.shape[0]):
                    for x in range(0, trans_div_matrix.shape[1]):
                        trans_pixel_value = trans_div_matrix[y, x]  # y,x are the row and column of the pixel in the individual trans matrix
                        abs_pixel_bin_x = chr_2_bin_abs_start + x
                        abs_pixel_bin_y = chr_1_bin_abs_start + y  # abs_pixel_bin_x, abs_pixel_bin_y are the coordinate of the pixel of the entire matrix
                        try:
                            chr_1_comp = next(tb.query(bin_dic[abs_pixel_bin_y][0], bin_dic[abs_pixel_bin_y][1], bin_dic[abs_pixel_bin_y][2]))[3]
                        except:
                            chr_1_comp = "chr_1_comp_nv"
                        try:
                            chr_2_comp = next(tb.query(bin_dic[abs_pixel_bin_x][0], bin_dic[abs_pixel_bin_x][1], bin_dic[abs_pixel_bin_x][2]))[3]
                        except:
                            chr_2_comp = "chr_2_comp_nv"
                        if chr_1_comp == chr_2_comp:
                            for trans_pixel_value_threshold in np.arange(min_trans_pixel_value, max_trans_pixel_value, increment_step):
                                if trans_pixel_value >= trans_pixel_value_threshold:
                                    csv_writer.writerow([bin_dic[abs_pixel_bin_y][0], bin_dic[abs_pixel_bin_x][0], chr_1_comp, trans_pixel_value_threshold, 'trans'])
                inter_time_end = time.perf_counter()
                print(f"trans matrix pixel between {chromosome_ls[i]} and {chromosome_ls[j]} finished, using {round(inter_time_end-inter_time_start, 2)} sec")
    print(f"{pixel_length_1} bp resolution pixel file containing subcompartment tagged value-thresholded pixels created")

    print(f"creating {increment_step} graded pixel subcompartment ratio file")
    sub_time_start = time.perf_counter()
    comp_tagged_df = pd.read_csv(div_matrix_pixel_file, delimiter='\t')
    with open(aggr_comp_pixel_count, 'w') as f5:
        csv_writer = csv.writer(f5, delimiter='\t')
        csv_writer.writerow(['compartment','interaction_type','counts','total_counts','ratio'])
        for cis_pixel_value_threshold in np.arange(min_cis_pixel_value, max_cis_pixel_value, increment_step):
            cis_comp_count_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0}
            for compartment in comp_ls:
                cis_comp_count_dic[compartment] = comp_tagged_df.loc[(comp_tagged_df['subcompartment'] == compartment) & (comp_tagged_df['interaction_type'] == 'cis') & (comp_tagged_df['threshold'] >= cis_pixel_value_threshold)].count()
                csv_writer.writerow([compartment, 'cis', cis_comp_count_dic[compartment], cis_pixel_total_counts_dic[compartment], round(100*cis_comp_count_dic[compartment]/cis_pixel_total_counts_dic[compartment], 2)])
        for trans_pixel_value_threshold in np.arange(min_trans_pixel_value, max_trans_pixel_value, increment_step):
            trans_comp_count_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0}
            for compartment in comp_ls:
                trans_comp_count_dic[compartment] = comp_tagged_df.loc[(comp_tagged_df['subcompartment'] == compartment) & (comp_tagged_df['interaction_type'] == 'trans') & (comp_tagged_df['threshold'] >= trans_pixel_value_threshold)].count()
                if compartment != 'B4':
                    csv_writer.writerow([compartment, 'trans', trans_comp_count_dic[compartment], trans_pixel_total_counts_dic[compartment], round(100*trans_comp_count_dic[compartment]/trans_pixel_total_counts_dic[compartment], 2)])
    sub_time_end = time.perf_counter()
    print(f"data aggregation finished in {round(sub_time_end - sub_time_start, 2)} seconds")
    end = time.perf_counter()
    print(f"Analysis finished in {round(end - start, 2)} seconds")


def main():
    parser = argparse.ArgumentParser(description= "To calculate pixel subcompartment contents of the division matrix of two cool files at certain resolution")
    parser.add_argument("-c1", help="cool file on top", dest="cool_1", type=str, required=True)
    parser.add_argument("-c2", help="cool file on bottom", dest="cool_2", type=str, required=True)
    parser.add_argument("-cf", help="output cis pixel value of the division matrix", dest="div_cis_pix", type=str, required=True)
    parser.add_argument("-tf", help="output trans pixel value of the division matrix", dest="div_trans_pix", type=str, required=True)
    parser.add_argument("-cp", help="subcompartment file", dest="comp_file", type=str, required=True)
    parser.add_argument("-is", help="increment step for the threshold", dest="incre_step", type=float, required=True)
    parser.add_argument("-ct", help="compartment tagged over-the-threshold pixel file", dest="div_matrix_thresholded_pix", type=str, required=True)
    parser.add_argument("-ap", help="aggregated over the threshold pixel compartment count ", dest="aggr_th_comp_ct", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()