#! /usr/bin/env python

import argparse
import cooler
import tabix
import numpy as np
import csv
import time

##  This program is used for calculating divsion matrix bin sum at the trans region
##  In addition, it will write pixels value of the third cooler file at the same location of interests.
##  polyA RNA expression information will also be added


def run(args):
    cool_file_1 = args.cool_1
    cool_file_2 = args.cool_2
    compartment_file = args.comp_file
    RNA_file = args.rna_file
    output_file = args.output
    comp_tb = tabix.open(compartment_file)
    polyA_RNA_tb = tabix.open(RNA_file)
    c1 = cooler.Cooler(cool_file_1)
    c2 = cooler.Cooler(cool_file_2)

    start = time.perf_counter()
    c1_matrix = c1.matrix(balance=False, sparse=False)[:]
    c2_matrix = c2.matrix(balance=False, sparse=False)[:]
    division_matrix = np.divide(c1_matrix, c2_matrix)
    division_matrix_noinf = np.where(np.isinf(division_matrix), np.nan, division_matrix)  # get rid of inf value
    division_matrix_noinf[division_matrix_noinf == 0] = np.nan # get rid of 0 value
    division_matrix_noinf_log = np.log10(division_matrix_noinf)

    chromosome_ls = []
    for i in range(1, 23):
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')

    # to get the bin and corresponding locations as a dictionary
    bins = c1.bins()[:]  # get the bin dataframe from one of the Hic matrix
    bin_dic = {}
    chrom_bin_dic = {}
    for chromosome in chromosome_ls:
        chrom_bins = bins.loc[bins['chrom'] == chromosome]
        chrom_bin_dic[chromosome] = (chrom_bins.index[0], chrom_bins.index[-1])   # chrom_bin_dic take the absolute bin for each chromosome
    for index, row in bins.iterrows():  # bins is a pandas dataframe containing the bin information from the cooler file
        bin_dic[index] = [row[0], row[1], row[2]]

    bin_max = int(chrom_bin_dic['chrX'][1])  # bin_max is the largest bin considered
    trans_bin_sum_dic = {}
    cis_bin_sum_dic = {}
    for bin_num in range(0, chrom_bin_dic['chrX'][1] + 1):  # initialize the bin_sum_dic to 0 for each bin
        trans_bin_sum_dic[bin_num] = 0
        cis_bin_sum_dic[bin_num] = 0

    for chromosome in chromosome_ls:  # without chrY and chrM
        chr_start_bin = chrom_bin_dic[chromosome][0]
        chr_end_bin = chrom_bin_dic[chromosome][1]
        print(f"chromosome {chromosome} 's bin start at {chr_start_bin} and end at {chr_end_bin}")
        for i in range(chr_start_bin, chr_end_bin + 1):
            for j in range(chr_start_bin, chr_end_bin + 1):
                if not np.isnan(division_matrix_noinf_log[i, j]):
                    cis_bin_sum_dic[i] += division_matrix_noinf_log[i, j]
            for j1 in range(0, chr_start_bin):
                if not np.isnan(division_matrix_noinf_log[i, j1]):
                    print(division_matrix_noinf_log[i, j1])
                    trans_bin_sum_dic[i] += division_matrix_noinf_log[i, j1]
            for j2 in range(chr_end_bin + 1, bin_max + 1):
                if not np.isnan(division_matrix_noinf_log[i, j2]):
                    trans_bin_sum_dic[i] += division_matrix_noinf_log[i, j2]

    with open(output_file, "w") as f1:
        csv_writer = csv.writer(f1, delimiter='\t')
        csv_writer.writerow(["bin_num", "chr_name", "chr_start", "chr_end", "subcompartment_start", "subcompartment_end", "subcompartment", "cis_bins_sum", "trans_bin_sum", "total_bin_sum", "gene_type", "fkpm", "gene_id"])  # write the header of the output file
        for bin_number in range(0, chrom_bin_dic['chrX'][1]+1):
            chr_name = bin_dic[bin_number][0]  ## chrom_bin_dic looks like chrX: 'chrX': (892, 939)
            chr_start = int(bin_dic[bin_number][1])
            chr_end = int(bin_dic[bin_number][2])
            cis_bin_sum = cis_bin_sum_dic[bin_number]
            trans_bin_sum = trans_bin_sum_dic[bin_number]
            bin_sum = cis_bin_sum + trans_bin_sum
            comp_query = comp_tb.query(chr_name, chr_start, chr_end)
            comp_query_ls = []
            for comp_query_result in comp_query:
                print(comp_query_result)
                comp_query_ls.append(comp_query_result)
            for comp_query_item in comp_query_ls:
                if int(comp_query_item[1]) <= chr_start:
                    comp_start = chr_start
                else:
                    comp_start = int(comp_query_item[1])

                if int(comp_query_item[2]) >= chr_end:
                    comp_end = chr_end
                else:
                    comp_end = int(comp_query_item[2])
                compartment = comp_query_item[3]
                polyA_RNA_query = polyA_RNA_tb.query(chr_name, comp_start, comp_end)
                for query_result in polyA_RNA_query:
                    if query_result[7] != "nan\r":  # row 8 of the query result contain FKPM information of the gene
                        gene_type = query_result[6].split(";")[2].split("=")[1]
                        gene_id = query_result[6].split(";")[1].split("=")[1]
                        fkpm_value = float(query_result[7])
                        csv_writer.writerow([bin_number, chr_name, chr_start, chr_end, comp_start, comp_end, compartment, cis_bin_sum, trans_bin_sum, bin_sum, gene_type, fkpm_value, gene_id])
    end = time.perf_counter()
    print(f"Analysis finished in {round(end - start, 2)} seconds")


def main():
    parser = argparse.ArgumentParser(description="To calculate pixel subcompartment contents of the division matrix of two cool files at certain resolution, and write the pixel value from the third cool file for comparison")
    parser.add_argument("-c1", help="cool file on top", dest="cool_1", type=str, required=True)
    parser.add_argument("-c2", help="cool file on bottom", dest="cool_2", type=str, required=True)
    parser.add_argument("-o", help="output tsv file containing the bin sum of the division matrix, together with", dest="output", type=str, required=True)
    parser.add_argument("-c", help="indexed subcompartment file", dest="comp_file", type=str, required=True)
    parser.add_argument("-r", help="indexed RNA expression file", dest="rna_file", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
