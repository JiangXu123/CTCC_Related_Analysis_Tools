#! /usr/bin/env python

import argparse
import cooler
import tabix
import numpy as np
import csv
import pandas as pd
import gzip
import time
from matplotlib import pyplot as plt
import seaborn as sns

'''This program is used for calculating the cool matrix pixel by pixel to find the pixel's contents in term of average subcompartment length
which is equal to the vertical pixel subcompartment length total plus horizontal pixel subcompartment length total, for each pixel, a table like beow
will be generated and write into the file:
chr_1_name  chr_1_start chr_1_end   chr_2_name  chr_2_start chr_2_end   A1  A1_average_length   division_value  third_matrix_value  cis
chr_1_name  chr_1_start chr_1_end   chr_2_name  chr_2_start chr_2_end   A2  A2_average_length   division_value  third_matrix_value  cis
chr_1_name  chr_1_start chr_1_end   chr_2_name  chr_2_start chr_2_end   B1  B1_average_length   division_value  third_matrix_value  cis
chr_1_name  chr_1_start chr_1_end   chr_2_name  chr_2_start chr_2_end   B2  B2_average_length   division_value  third_matrix_value  cis
chr_1_name  chr_1_start chr_1_end   chr_2_name  chr_2_start chr_2_end   B3  B3_average_length   division_value  third_matrix_value  cis
enrichment = sum_of_pixels_covered_average_subcompartment_length / total_subcompartment_length'''


def run(args):
    cool_file_1 = args.cool_1
    compartment_file = args.comp_file
    matrix_pixel_comp_file = args.mat_comp_pix
    aggr_comp_pixel_file = args.aggr_th_comp_ct
    comp_length_pixel_threshold_fig_file = args.fig

    pixel_length = int(cool_file_1.split("/")[-1])
    tb = tabix.open(compartment_file)
    c1 = cooler.Cooler(cool_file_1)
    start = time.perf_counter()
    # generate a chromosome list to contain all the chromosomes
    chromosome_ls = []
    for i in range(1, 23):
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')
    chromosome_ls.append('chrY')
    chromosome_ls.append('chrM')
    # the chromosome list is to contain all the chromosomes that is listed in the compartment file

    # the first step is to calculate how many pixels are same compartment pixel for a certain subcompartment
    # comp_len_dic contains the length of each subcompartment (in bp) of the GM12878 cell line
    # cis_pixel_total_counts_dic contains the theoretical cis same compartment pixel counts for each subcompartment,
    # used to normalize the data considering the difference in compartments' length
    comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
    cis_comp_ls = ["A1", "A2", "B1", "B2", "B3", "B4", "NA"]
    trans_comp_ls = ["A1", "A2", "B1", "B2", "B3", "NA"]
    with gzip.open(compartment_file, 'rt') as f1:
        csv_reader = csv.reader(f1, delimiter='\t')
        for line in csv_reader:
            comp_len_dic[line[3]] += (int(line[2]) - int(line[1]))  # generate each compartment length total
    print(comp_len_dic)
    # the last step is to parse the cis and trans matrix to get the pixels and their belonged compartments,
    # only same compartment pixels are considered.

    # to get the bin and corresponding locations as a dictionary
    bins = c1.bins()[:]  # get the bin dataframe from one of the Hic matrix
    bin_dic = {}
    for index, row in bins.iterrows():  # bins is a pandas dataframe containing the bin information from the cooler file
        bin_dic[index] = [row[0], row[1], row[2]]
    print(f"bins for the matrix at resolution {pixel_length} have been generated")

    chrom_exclusion_ls = ['chrM', 'chrY']  # chrM has significant higher pixel values, that is excluded from the analysis

    with open(matrix_pixel_comp_file, "w") as f2:
        csv_writer = csv.writer(f2, delimiter='\t')
        csv_writer.writerow(['chr_1_name', 'chr_1_start', 'chr_1_end', 'chr_2_name', 'chr_2_start', 'chr_2_end', 'subcompartment', 'subcompartment_length', 'pixel_value', 'interaction_type'])  # write the header of the output file
        # 'subcompartment_length' refers to pixels_covered_average_subcompartment_length
        for k in range(0, len(chromosome_ls)):
            inter_time_start = time.perf_counter()
            chr_name = chromosome_ls[k]
            if not (chr_name in chrom_exclusion_ls):
                cis_1 = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[k])  # cis_1, cis_2, and cis_3 are the cis matrix of the corresponding chromosome of three cool files, respectively.
                chr_bins = bins.loc[bins['chrom'] == chr_name]
                chr_bin_abs_start = chr_bins.index[0]  # to get the absolute start bin position of the chromosome
                for r in range(0, cis_1.shape[0]):
                    for c in range(0, cis_1.shape[1]):
                        cis_pixel_value = cis_1[r, c]  # r, c are the row and column of the pixel in the individual trans matrix
                        r_pixel_comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
                        c_pixel_comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
                        average_pixel_comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
                        if (not np.isinf(cis_pixel_value)) & (not np.isnan(cis_pixel_value)) & (cis_pixel_value != 0):
                            abs_pixel_bin_c = chr_bin_abs_start + c
                            abs_pixel_bin_r = chr_bin_abs_start + r

                            cis_chr_1_name = bin_dic[abs_pixel_bin_r][0]
                            cis_chr_1_start = bin_dic[abs_pixel_bin_r][1]
                            cis_chr_1_end = bin_dic[abs_pixel_bin_r][2]

                            cis_chr_2_name = bin_dic[abs_pixel_bin_c][0]
                            cis_chr_2_start = bin_dic[abs_pixel_bin_c][1]
                            cis_chr_2_end = bin_dic[abs_pixel_bin_c][2]

                            chr_1_query_ls = []
                            chr_2_query_ls = []

                            try:
                                chr_1_query_results = tb.query(cis_chr_1_name, cis_chr_1_start, cis_chr_1_end)
                                for result_1 in chr_1_query_results:
                                    chr_1_query_ls.append(result_1)  # For low resolution pixel (such as 5 M/pixel), a pixel may encompass several subcompartment,
                            except:
                                pass

                            try:
                                chr_2_query_results = tb.query(cis_chr_2_name, cis_chr_2_start, cis_chr_2_end)
                                for result_2 in chr_2_query_results:
                                    chr_2_query_ls.append(result_2)
                            except:
                                pass

                            for comp_q_1 in chr_1_query_ls:  # count the pixel contained subcompartment length
                                if int(comp_q_1[1]) <= cis_chr_1_start:
                                    comp_1_start = cis_chr_1_start
                                else:
                                    comp_1_start = int(comp_q_1[1])

                                if int(comp_q_1[2]) >= cis_chr_1_end:
                                    comp_1_end = cis_chr_1_end
                                else:
                                    comp_1_end = int(comp_q_1[2])

                                comp_1_length = comp_1_end - comp_1_start
                                comp_1_compartment = comp_q_1[3]
                                r_pixel_comp_len_dic[comp_1_compartment] += comp_1_length

                            for comp_q_2 in chr_2_query_ls:
                                if int(comp_q_2[1]) <= cis_chr_2_start:
                                    comp_2_start = cis_chr_2_start
                                else:
                                    comp_2_start = int(comp_q_2[1])

                                if int(comp_q_2[2]) >= cis_chr_2_end:
                                    comp_2_end = cis_chr_2_end
                                else:
                                    comp_2_end = int(comp_q_2[2])

                                comp_2_length = comp_2_end - comp_2_start
                                comp_2_compartment = comp_q_2[3]
                                c_pixel_comp_len_dic[comp_2_compartment] += comp_2_length

                            for comp in cis_comp_ls:
                                average_pixel_comp_len_dic[comp] = (r_pixel_comp_len_dic[comp] + c_pixel_comp_len_dic[comp]) / 2
                                csv_writer.writerow([cis_chr_1_name, cis_chr_1_start, cis_chr_1_end, cis_chr_2_name, cis_chr_2_start, cis_chr_2_end, comp, average_pixel_comp_len_dic[comp], cis_pixel_value, 'cis'])
                inter_time_end = time.perf_counter()
            print(f"{chromosome_ls[k]} cis matrix parsed in {round(inter_time_end - inter_time_start, 2)} seconds")

        for i in range(0, len(chromosome_ls)):  # generate the over-the-threshold trans pixel contents first
            for j in range(i + 1, len(chromosome_ls)):
                inter_time_start = time.perf_counter()
                chr_name_1 = chromosome_ls[i]
                chr_name_2 = chromosome_ls[j]
                if (not (chr_name_1 in chrom_exclusion_ls)) & (not (chr_name_2 in chrom_exclusion_ls)):
                    trans_1 = c1.matrix(sparse=False, balance=True).fetch(chromosome_ls[i], chromosome_ls[j])
                    chr_bins_1 = bins.loc[bins['chrom'] == chr_name_1]
                    chr_bins_2 = bins.loc[bins['chrom'] == chr_name_2]
                    chr_1_bin_abs_start = chr_bins_1.index[0]  # get the relative starting bin number for the chromosome
                    chr_2_bin_abs_start = chr_bins_2.index[0]
                    for y in range(0, trans_1.shape[0]):
                        for x in range(0, trans_1.shape[1]):
                            trans_pixel_value = trans_1[y, x]  # y, x are the row and column of the pixel in the individual trans matrix
                            r_pixel_comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
                            c_pixel_comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
                            average_pixel_comp_len_dic = {"A1": 0, "A2": 0, "B1": 0, "B2": 0, "B3": 0, "B4": 0, "NA": 0}
                            if (not np.isinf(trans_pixel_value)) & (not np.isnan(trans_pixel_value)) & (trans_pixel_value != 0):
                                abs_pixel_bin_x = chr_2_bin_abs_start + x
                                abs_pixel_bin_y = chr_1_bin_abs_start + y  # abs_pixel_bin_x, abs_pixel_bin_y are the coordinate of the pixel of the entire matrix
                                trans_chr_1_name = bin_dic[abs_pixel_bin_y][0]
                                trans_chr_1_start = bin_dic[abs_pixel_bin_y][1]
                                trans_chr_1_end = bin_dic[abs_pixel_bin_y][2]
                                trans_chr_2_name = bin_dic[abs_pixel_bin_x][0]
                                trans_chr_2_start = bin_dic[abs_pixel_bin_x][1]
                                trans_chr_2_end = bin_dic[abs_pixel_bin_x][2]
                                chr_1_query_ls = []
                                chr_2_query_ls = []

                                try:
                                    chr_1_query_results = tb.query(trans_chr_1_name, trans_chr_1_start, trans_chr_1_end)
                                    for result_1 in chr_1_query_results:
                                        chr_1_query_ls.append(result_1)
                                except:
                                    pass
                                try:
                                    chr_2_query_results = tb.query(trans_chr_2_name, trans_chr_2_start, trans_chr_2_end)
                                    for result_2 in chr_2_query_results:
                                        chr_2_query_ls.append(result_2)
                                except:
                                    pass

                                for comp_q_1 in chr_1_query_ls:  # if both query results contains the same compartment, then it's regarded as a effective same compartment contact
                                    if int(comp_q_1[1]) <= trans_chr_1_start:
                                        comp_1_start = trans_chr_1_start
                                    else:
                                        comp_1_start = int(comp_q_1[1])
                                    if int(comp_q_1[2]) >= trans_chr_1_end:
                                        comp_1_end = trans_chr_1_end
                                    else:
                                        comp_1_end = int(comp_q_1[2])
                                    comp_1_length = comp_1_end - comp_1_start
                                    comp_1_compartment = comp_q_1[3]
                                    r_pixel_comp_len_dic[comp_1_compartment] += comp_1_length

                                for comp_q_2 in chr_2_query_ls:
                                    if int(comp_q_2[1]) <= trans_chr_2_start:
                                        comp_2_start = trans_chr_2_start
                                    else:
                                        comp_2_start = int(comp_q_2[1])
                                    if int(comp_q_2[2]) >= trans_chr_2_end:
                                        comp_2_end = trans_chr_2_end
                                    else:
                                        comp_2_end = int(comp_q_2[2])
                                    comp_2_length = comp_2_end - comp_2_start
                                    comp_2_compartment = comp_q_2[3]
                                    c_pixel_comp_len_dic[comp_2_compartment] += comp_2_length

                                for comp in trans_comp_ls:
                                    average_pixel_comp_len_dic[comp] = (r_pixel_comp_len_dic[comp] + c_pixel_comp_len_dic[comp]) / 2
                                    csv_writer.writerow([trans_chr_1_name, trans_chr_1_start, trans_chr_1_end, trans_chr_2_name, trans_chr_2_start, trans_chr_2_end, comp, average_pixel_comp_len_dic[comp], trans_pixel_value, 'trans'])
                inter_time_end = time.perf_counter()
                print(f"trans matrix pixel between {chromosome_ls[i]} and {chromosome_ls[j]} finished, using {round(inter_time_end - inter_time_start, 2)} sec")
    print(f"{pixel_length} bp resolution pixel file containing subcompartment tagged value-thresholded pixels created")

    # the second step is to calculate the maximum value and the minimum value of the filtered pixel
    pixel_df = pd.read_csv(matrix_pixel_comp_file, delimiter='\t')
    cis_pixel_min = pixel_df.loc[pixel_df['interaction_type'] == 'cis'].describe().loc['min']['pixel_value']
    cis_pixel_max = pixel_df.loc[pixel_df['interaction_type'] == 'cis'].describe().loc['max']['pixel_value']
    trans_pixel_min = pixel_df.loc[pixel_df['interaction_type'] == 'trans'].describe().loc['min']['pixel_value']
    trans_pixel_max = pixel_df.loc[pixel_df['interaction_type'] == 'trans'].describe().loc['max']['pixel_value']

    print("All values are calculated without chrM and chrY")
    print(f"The highest non-inf cis value of the contact matrix is {cis_pixel_max}")
    print(f"The lowest non-zero cis value of the contact matrix is {cis_pixel_min}")
    print(f"The highest non-inf trans value of the division matrix is {trans_pixel_max}")
    print(f"The lowest non-zero trans value of the division matrix is {trans_pixel_min}")

    # the last step is to create the aggregated counts for filtered pixel counts for each subcompartment
    print(f"creating 0.1 spaced thresholded pixel value with subcompartment ratio file")
    sub_time_start = time.perf_counter()

    with open(aggr_comp_pixel_file, 'w') as f3:
        csv_writer = csv.writer(f3, delimiter='\t')
        csv_writer.writerow(['compartment', 'interaction_type', 'threshold', 'length', 'total_length', 'ratio'])
        for cis_pixel_value_threshold in np.linspace(0, cis_pixel_max, 100):  # the value of each item in the division matrix is somewhere between 0 and positive infinity
            cis_comp_length_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0, 'B4': 0}
            for compartment in cis_comp_ls:
                cis_comp_length_dic[compartment] = pixel_df.loc[(pixel_df['subcompartment'] == compartment) & (pixel_df['interaction_type'] == 'cis') & (pixel_df['pixel_value'] >= cis_pixel_value_threshold)]['subcompartment_length'].sum()
                csv_writer.writerow([compartment, 'cis', cis_pixel_value_threshold, cis_comp_length_dic[compartment], comp_len_dic[compartment], round(100 * cis_comp_length_dic[compartment] / comp_len_dic[compartment], 5)])
        for trans_pixel_value_threshold in np.linspace(0, trans_pixel_max, 100):
            trans_comp_length_dic = {'A1': 0, 'A2': 0, 'B1': 0, 'B2': 0, 'B3': 0}
            for compartment in trans_comp_ls:
                trans_comp_length_dic[compartment] = pixel_df.loc[(pixel_df['subcompartment'] == compartment) & (pixel_df['interaction_type'] == 'trans') & (pixel_df['pixel_value'] >= trans_pixel_value_threshold)]['subcompartment_length'].sum()
                csv_writer.writerow([compartment, 'trans', trans_pixel_value_threshold, trans_comp_length_dic[compartment], comp_len_dic[compartment], round(100 * trans_comp_length_dic[compartment] / comp_len_dic[compartment], 5)])
    sub_time_end = time.perf_counter()
    print(f"data aggregation finished in {round(sub_time_end - sub_time_start, 2)} seconds")
    end = time.perf_counter()
    print(f"Analysis finished in {round(end - start, 2)} seconds")
    # cis_comp_length_dic contains the sum of averaged pixel compartment length exceeding graduated values (from 0 to div_cis_pixel_max)
    # trans_comp_length_dic contains the sum of averaged pixel compartment length exceeding graduated values (from 0 to div_cis_pixel_max)

    # plot the aggregated file in a 2 rows by 2 column fashion.
    # B4 subcompartment was not included in the plot.
    aggr_pixel_df = pd.read_csv(aggr_comp_pixel_file, delimiter='\t')
    cis_enrichment_fold_max = aggr_pixel_df.loc[(aggr_pixel_df['interaction_type'] == "cis") & (aggr_pixel_df['compartment'] != 'B4')].describe().loc['max']['ratio']
    trans_enrichment_fold_max = aggr_pixel_df.loc[(aggr_pixel_df['interaction_type'] == "trans") & (aggr_pixel_df['compartment'] != 'B4')].describe().loc['max']['ratio']

    cis_df = aggr_pixel_df.loc[(aggr_pixel_df['interaction_type'] == "cis") & (aggr_pixel_df['compartment'] != 'B4')]
    trans_df = aggr_pixel_df.loc[(aggr_pixel_df['interaction_type'] == "trans") & (aggr_pixel_df['compartment'] != 'B4')]

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 8), dpi=300)  # size of the plot is 16 inch wide and 8 inch height, encomparsing two picture horizontally
    fig.tight_layout(pad=5)  # make enough space between subplots
    sns.lineplot(x='threshold', y='ratio', hue='compartment', palette=["#D55E00", "#CC79A7", "#999999", "#56B4E9", "#F0E442"], linewidth=3, data=cis_df, ax=ax[0])
    ax[0].set_yscale('log')
    ax[0].set_xlabel('Pixel Threshold', fontsize=20)
    ax[0].set_ylabel('Enrichment Fold', fontsize=20)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[0].text((cis_pixel_max/2), cis_enrichment_fold_max, "Cis", fontsize=20)
    ax[0].legend(frameon=False, fontsize=20)
    ax[0].tick_params(labelsize=20)

    sns.lineplot(x='threshold', y='ratio', hue='compartment', palette=["#D55E00", "#CC79A7", "#999999", "#56B4E9", "#F0E442"], linewidth=3, data=trans_df, ax=ax[1])
    ax[1].set_yscale('log')
    ax[1].set_xlabel('Pixel Threshold', fontsize=20)
    ax[1].set_ylabel('Enrichment Fold', fontsize=20)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].text((trans_pixel_max/2), trans_enrichment_fold_max, "Trans", fontsize=20)
    ax[1].legend(frameon=False, fontsize=20)
    ax[1].tick_params(labelsize=20)

    fig.savefig(comp_length_pixel_threshold_fig_file, format='png', dpi=300)

def main():
    parser = argparse.ArgumentParser(description="To calculate pixel subcompartment contents of the cool files at certain resolution and plot the enrichment fold-pixel_threshold figure")
    parser.add_argument("-co", help="cool file", dest="cool_1", type=str, required=True)
    parser.add_argument("-cp", help="subcompartment file", dest="comp_file", type=str, required=True)
    parser.add_argument("-pc", help="output file containing graded pixel value and subcompartments", dest="mat_comp_pix", type=str, required=True)
    parser.add_argument("-ap", help="aggregated over the threshold pixel compartment count ", dest="aggr_th_comp_ct", type=str, required=True)
    parser.add_argument("-f", help="output lineplot figure file, should be ended with .png", dest="fig", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
