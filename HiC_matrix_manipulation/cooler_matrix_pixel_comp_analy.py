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
    cool_file = args.cool
    pixel_value_comp_file = args.pix_comp
    compartment_file = args.comp_file
    aggr_comp_pixel_freq = args.aggr_th_comp_ct

    pixel_length = int(cool_file.split("/")[-1])
    tb = tabix.open(compartment_file)
    c1 = cooler.Cooler(cool_file)
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

    # comp_len_dic contains the length of subcompartments (bp) of the GM12878 cell line
    # cis_pixel_total_counts_dic contains the cis same compartment pixel counts for each subcompartment,
    # used to normalize the data considering the difference in compartments' length
    comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}

    with gzip.open(compartment_file, 'rt') as f1:
        csv_reader = csv.reader(f1, delimiter='\t')
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
        print(chr_comp_dic[chromosome]['A1'])
        cis_pixel_total_counts_dic['A1'] += int(chr_comp_dic[chromosome]['A1'] / pixel_length) * int(chr_comp_dic[chromosome]['A1'] / pixel_length)
        cis_pixel_total_counts_dic['A2'] += int(chr_comp_dic[chromosome]['A2'] / pixel_length) * int(chr_comp_dic[chromosome]['A2'] / pixel_length)
        cis_pixel_total_counts_dic['B1'] += int(chr_comp_dic[chromosome]['B1'] / pixel_length) * int(chr_comp_dic[chromosome]['B1'] / pixel_length)
        cis_pixel_total_counts_dic['B2'] += int(chr_comp_dic[chromosome]['B2'] / pixel_length) * int(chr_comp_dic[chromosome]['B2'] / pixel_length)
        cis_pixel_total_counts_dic['B3'] += int(chr_comp_dic[chromosome]['B3'] / pixel_length) * int(chr_comp_dic[chromosome]['B3'] / pixel_length)
        cis_pixel_total_counts_dic['B4'] += int(chr_comp_dic[chromosome]['B4'] / pixel_length) * int(chr_comp_dic[chromosome]['B4'] / pixel_length)
        cis_pixel_total_counts_dic['NA'] += int(chr_comp_dic[chromosome]['NA'] / pixel_length) * int(chr_comp_dic[chromosome]['NA'] / pixel_length)
    print(f'the reference cis_pixel_total_same_compartment_pixel dictionary is {cis_pixel_total_counts_dic}')

    total_comp_pixel_number_dic = {}  # total_comp_pixel_number_dic = trans_pixel_total_counts_dic + cis_pixel_total_counts_dic
    trans_pixel_total_counts_dic = {}
    comp_ls = ['A1', 'A2', 'B1', 'B2', 'B3', 'B4', 'NA']

    # to calculate the total pixel number belong to the same compartment interaction matrix, which include trans and cis
    for compartment in comp_ls:
        total_comp_pixel_number_dic[compartment] = int(comp_len_dic[compartment] / pixel_length) * int(comp_len_dic[compartment] / pixel_length) # this is a rough estimation of the subcompartment-involved pixel counts
        trans_pixel_total_counts_dic[compartment] = total_comp_pixel_number_dic[compartment] - cis_pixel_total_counts_dic[compartment]
    print(f'the reference total_pixel_total_same_compartment_pixel dictionary is {total_comp_pixel_number_dic}')
    print(f'the reference trans_pixel_total_same_compartment_pixel dictionary is {trans_pixel_total_counts_dic}')
    # the last step is to parse the cis and trans matrix to get the pixels and their belonged compartments,
    # only same compartment pixels are considered.

    # to get the bin and corresponding locations as a dictionary
    bins = c1.bins()[:]  # get the bin dataframe from one of the Hic matrix
    bin_dic = {}
    for index, row in bins.iterrows():  # bins is a pandas dataframe containing the bin information from the cooler file
        bin_dic[index] = [row[0], row[1], row[2]]
    print(f"bins for the matrix at resolution {pixel_length} have been generated")

    chrom_exclusion_ls = ['chrM', 'chrY']
    with open(pixel_value_comp_file, "w") as f2:
        csv_writer = csv.writer(f2, delimiter='\t')
        csv_writer.writerow(['chr_1_name', 'chr_1_start', 'chr_1_end', 'chr_2_name', 'chr_2_start', 'chr_2_end', 'subcompartment', 'pixel_value', 'interaction_type'])  # write the header of the output file
        for k in range(0, len(chromosome_ls)):
            inter_time_start = time.perf_counter()
            cis = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[k])
            chr_name = chromosome_ls[k]  # chromosome name
            chr_bins = bins.loc[bins['chrom'] == chr_name]
            chr_bin_abs_start = chr_bins.index[0]
            for r in range(0, cis.shape[0]):
                for c in range(0, cis.shape[1]):
                    abs_pixel_bin_c = chr_bin_abs_start + c
                    abs_pixel_bin_r = chr_bin_abs_start + r
                    if not (chr_name in chrom_exclusion_ls): # chrM has significant higher pixel values, that is excluded from the analysis
                        cis_pixel_value = cis[r, c]  # r, c are the row and column of the pixel in the individual trans matrix
                        if (not np.isinf(cis_pixel_value)) & (not np.isnan(cis_pixel_value)) & (cis_pixel_value != 0):
                            chr_1_start = bin_dic[abs_pixel_bin_r][1]
                            chr_1_end = bin_dic[abs_pixel_bin_r][2]
                            chr_2_start = bin_dic[abs_pixel_bin_r][1]
                            chr_2_end = bin_dic[abs_pixel_bin_r][2]
                            try:
                                r_comp = next(tb.query(bin_dic[abs_pixel_bin_r][0], bin_dic[abs_pixel_bin_r][1], bin_dic[abs_pixel_bin_r][2]))[3]
                            except:
                                r_comp = "r_comp_nv"
                            try:
                                c_comp = next(tb.query(bin_dic[abs_pixel_bin_c][0], bin_dic[abs_pixel_bin_c][1], bin_dic[abs_pixel_bin_c][2]))[3]
                            except:
                                c_comp = "c_comp_nv"
                            if r_comp == c_comp:
                                csv_writer.writerow([chr_name, chr_1_start, chr_1_end, chr_name, chr_2_start, chr_2_end, r_comp, cis_pixel_value, 'cis'])
            inter_time_end = time.perf_counter()
            print(f"{chromosome_ls[k]} cis matrix parsed in {round(inter_time_end - inter_time_start, 2)} seconds")
        for i in range(0, len(chromosome_ls)):  # generate the over-the-threshold trans pixel contents first
            for j in range(i + 1, len(chromosome_ls)):
                inter_time_start = time.perf_counter()
                chr_name_1 = chromosome_ls[i]  # first chromsome name
                chr_name_2 = chromosome_ls[j]  # second chromosome name
                if (not (chr_name_1 in chrom_exclusion_ls)) & (not (chr_name_2 in chrom_exclusion_ls)):  # chrM with other chromosomes have significant higher pixel values, that should be excluded from the analysis
                    trans = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[i], chromosome_ls[j])
                    chr_bins_1 = bins.loc[bins['chrom'] == chr_name_1]
                    chr_bins_2 = bins.loc[bins['chrom'] == chr_name_2]
                    chr_1_bin_abs_start = chr_bins_1.index[0]
                    chr_2_bin_abs_start = chr_bins_2.index[0]
                    for y in range(0, trans.shape[0]):
                        for x in range(0, trans.shape[1]):
                            trans_pixel_value = trans[y, x]  # y,x are the row and column of the pixel in the individual trans matrix
                            if (not np.isinf(trans_pixel_value)) & (not np.isnan(trans_pixel_value)) & (trans_pixel_value != 0):
                                abs_pixel_bin_x = chr_2_bin_abs_start + x
                                abs_pixel_bin_y = chr_1_bin_abs_start + y  # abs_pixel_bin_x, abs_pixel_bin_y are the coordinate of the pixel of the entire matrix
                                trans_chr_1_start = bin_dic[abs_pixel_bin_y][1]
                                trans_chr_1_end = bin_dic[abs_pixel_bin_y][2]
                                trans_chr_2_start = bin_dic[abs_pixel_bin_x][1]
                                trans_chr_2_end = bin_dic[abs_pixel_bin_x][2]
                                try:
                                    chr_1_comp = next(tb.query(bin_dic[abs_pixel_bin_y][0], bin_dic[abs_pixel_bin_y][1], bin_dic[abs_pixel_bin_y][2]))[3]
                                except:
                                    chr_1_comp = "chr_1_comp_nv"
                                try:
                                    chr_2_comp = next(tb.query(bin_dic[abs_pixel_bin_x][0], bin_dic[abs_pixel_bin_x][1], bin_dic[abs_pixel_bin_x][2]))[3]
                                except:
                                    chr_2_comp = "chr_2_comp_nv"
                                if chr_1_comp == chr_2_comp:
                                    csv_writer.writerow([chr_name_1, trans_chr_1_start, trans_chr_1_end, chr_name_2, trans_chr_2_start, trans_chr_2_end, chr_1_comp, trans_pixel_value, 'trans'])
                inter_time_end = time.perf_counter()
                print(f"trans matrix pixel between {chromosome_ls[i]} and {chromosome_ls[j]} finished, using {round(inter_time_end - inter_time_start, 2)} sec")
    print(f"{pixel_length} bp resolution pixel file containing subcompartment tagged pixel value created")

    ### the second step is to calculate the maximum value and the minimum value of the filtered pixel
    pixel_df = pd.read_csv(pixel_value_comp_file, delimiter='\t')
    cis_pixel_max = pixel_df.loc[pixel_df['interaction_type'] == 'cis'].describe().loc['max']['pixel_value']
    trans_pixel_max = pixel_df.loc[pixel_df['interaction_type'] == 'trans'].describe().loc['max']['pixel_value']
    print(f"The highest non-inf cis value of the division matrix is {cis_pixel_max}")
    print(f"The highest non-inf trans value of the division matrix is {trans_pixel_max}")

    ### the last step is to creat the aggregated counts for filtered pixel counts for each subcompartment
    sub_time_start = time.perf_counter()

    with open(aggr_comp_pixel_freq, 'w') as f3:
        csv_writer = csv.writer(f3, delimiter='\t')
        csv_writer.writerow(['compartment', 'interaction_type', 'threshold', 'counts', 'total_counts', 'ratio'])
        for cis_pixel_value_threshold in np.linspace(0, cis_pixel_max, 100):
            cis_comp_count_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0}
            for compartment in comp_ls:
                cis_comp_count_dic[compartment] = pixel_df.loc[(pixel_df['subcompartment'] == compartment) & (pixel_df['interaction_type'] == 'cis') & (pixel_df['pixel_value'] >= cis_pixel_value_threshold)].count()['interaction_type']
                print(f'cis pixels over {cis_pixel_value_threshold} and has same {compartment} has  {cis_comp_count_dic[compartment]} count')
                csv_writer.writerow([compartment, 'cis', cis_pixel_value_threshold, cis_comp_count_dic[compartment], cis_pixel_total_counts_dic[compartment], round(100 * cis_comp_count_dic[compartment] / cis_pixel_total_counts_dic[compartment], 2)])
        for trans_pixel_value_threshold in np.linspace(0, trans_pixel_max, 100):
            trans_comp_count_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0}
            for compartment in comp_ls:
                trans_comp_count_dic[compartment] = pixel_df.loc[(pixel_df['subcompartment'] == compartment) & (pixel_df['interaction_type'] == 'trans') & (pixel_df['pixel_value'] >= trans_pixel_value_threshold)].count()['interaction_type']
                print(f'trans pixels over {trans_pixel_value_threshold} and belongs to same {compartment} has  {trans_comp_count_dic[compartment]} counts')
                if compartment != 'B4':
                    csv_writer.writerow([compartment, 'trans', trans_pixel_value_threshold, trans_comp_count_dic[compartment], trans_pixel_total_counts_dic[compartment], round(100 * trans_comp_count_dic[compartment] / trans_pixel_total_counts_dic[compartment], 6)])
    sub_time_end = time.perf_counter()
    print(f"data aggregation finished in {round(sub_time_end - sub_time_start, 2)} seconds")
    end = time.perf_counter()
    print(f"Analysis finished in {round(end - start, 2)} seconds")


def main():
    parser = argparse.ArgumentParser(description="To calculate pixel subcompartment contents of the hic matrix of two cool files at certain resolution")
    parser.add_argument("-c1", help="cool file to be analyzed", dest="cool", type=str, required=True)
    parser.add_argument("-cp", help="subcompartment file", dest="comp_file", type=str, required=True)
    parser.add_argument("-pc", help="output file containing graded pixel value and subcompartment information", dest="pix_comp", type=str, required=True)
    parser.add_argument("-ap", help="aggregated over the threshold pixel compartment count", dest="aggr_th_comp_ct", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()