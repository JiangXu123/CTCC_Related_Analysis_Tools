#! /usr/bin/env python

import csv
import time
from operator import itemgetter
import argparse
import ast
import pandas as pd


def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file_prefix = args.out_pf
    same_bc_reads_count_lower_thresh_hold = args.sbrclt
    same_bc_reads_count_higher_thresh_hold = args.sbrcht
    smaller_output_file = output_file_prefix + "sbrlt_" + str(same_bc_reads_count_lower_thresh_hold) + "_l.pair"
    larger_output_file = output_file_prefix + "sbrlt_" + str(same_bc_reads_count_lower_thresh_hold) + "_sbrht_" + str(same_bc_reads_count_higher_thresh_hold) + ".pair"
    same_bc_read_list_file = output_file_prefix + "_sbcr.tsv"
    same_bc_read_length_list_file = output_file_prefix + "_sbcrl.tsv"
    mapping_thresh_hold = args.mpqt
    output_pair_larger_than_same_bc_read_threshold = args.lp
    data_ls = []
    # get same barcoded reads to a list and sort it according to the barcode
    if not mapping_thresh_hold:   # if no mapping quality threshold was inputted
        mapping_thresh_hold = 0   # default mapping quality threshold is 0
    if not output_pair_larger_than_same_bc_read_threshold:  # if output_pair_larger_than_same_bc_read_threshold option was not inputted as 'yes', do not output such pairs file,
        output_pair_larger_than_same_bc_read_threshold = 'no'
    with open(input_file, 'r') as csv_file1:
        csv_reader1 = csv.reader(csv_file1, delimiter='\t')
        for line in csv_reader1:
            if (int(line[4]) > 101010101) & (line[4][-2:] != '00') & (line[4][-4:-2] != '00') & (line[4][-6:-4] != '00') & (line[4][-8:-6] != '00') & (line[4][0:-8] != '00') & (line[1] != '') & (int(line[3]) >= mapping_thresh_hold):
                data_ls.append([line[1], line[2], line[4]])  # line[0] is read name, line[1] is mapping chromsome name, line[2] is position, line[4] is barcode
    sorted_ls = sorted(data_ls, key=itemgetter(2))  # sort the list based on the barcode, so same same barcoded reads will be following each other.
    qualified_read_count = len(sorted_ls)

    print('mapped_bc_reads sorting finished')
    ls = []
    all_ls = []
    temp = sorted_ls[0]

    ls.append(temp)     # code to cluster the reads with same barcode to a list 'ls'
    for i in range(1, len(sorted_ls)):   # and append the list to another list 'all_ls'
        if (temp[2] == sorted_ls[i][2]) & (sorted_ls[i] not in ls): # if neighboring two reads share same barcode and are not PCR duplicate
            ls.append(sorted_ls[i])
        else:
            if len(ls) > 1:  # if neighboring two reads don't share the same barcode and if the ls holding more than 1 same barcoded reads(2 or more)
                all_ls.append(ls)
            temp = sorted_ls[i]
            ls = []
            ls.append(temp)
    print('mapped_bc_with_same_barcode list generated')

    with open(same_bc_read_list_file, 'w') as csv_file2:
        csv_writer1 = csv.writer(csv_file2, delimiter='\t')
        csv_writer1.writerows(all_ls)
    print('same barcoded reads list generated')

    item_count = 0
    same_bc_reads_length_ls = []
    with open(same_bc_read_list_file, 'r') as csv_file3:
        csv_reader3 = csv.reader(csv_file3, delimiter='\t')
        with open(smaller_output_file, 'w') as csv_file4:
            csv_writer2 = csv.writer(csv_file4, delimiter='\t')
            with open(larger_output_file, 'w') as csv_file5:
                csv_writer3 = csv.writer(csv_file5, delimiter='\t')
                for item in csv_reader3:
                    process_start = time.perf_counter()
                    smaller_temp_ls = []
                    larger_temp_ls = []
                    same_bc_read_length = len(item)
                    same_bc_reads_length_ls.append(same_bc_read_length)
                    if same_bc_read_length <= same_bc_reads_count_lower_thresh_hold:
                        for indiv_ls in item:
                            smaller_temp_ls.append(ast.literal_eval(indiv_ls)) # convert each item in the list from a string to a list data format
                        for a in range(0, len(smaller_temp_ls) - 1):
                            for b in range(a + 1, len(smaller_temp_ls)):
                                csv_writer2.writerow([smaller_temp_ls[a][0], smaller_temp_ls[a][1], smaller_temp_ls[b][0], smaller_temp_ls[b][1], 1])
                    if output_pair_larger_than_same_bc_read_threshold == 'yes':
                        if (same_bc_read_length > same_bc_reads_count_lower_thresh_hold) & (same_bc_read_length <= same_bc_reads_count_higher_thresh_hold):
                            for indiv_ls in item:
                                larger_temp_ls.append(ast.literal_eval(indiv_ls))  # convert each item in the list from a string to a list data format
                            for c in range(0, len(larger_temp_ls) - 1):
                                for d in range(c + 1, len(larger_temp_ls)):
                                    csv_writer3.writerow([larger_temp_ls[c][0], larger_temp_ls[c][1], larger_temp_ls[d][0], larger_temp_ls[d][1], 1])
                    item_count += 1
                    process_end = time.perf_counter()
                    print(f'{item_count} cluster with same barcoded reads processed')
                    print(f'{item_count} cluster takes {round(process_end-process_start, 2)} seconds')
    same_bc_read_count_df = pd.DataFrame(same_bc_reads_length_ls, columns=['same_bc_read_number'])
    same_bc_read_count_df['count'] = 1
    aggr_same_bc_read_count_df = same_bc_read_count_df.groupby(['same_bc_read_number']).count()['count'].reset_index()
    aggr_same_bc_read_count_df.to_csv(same_bc_read_length_list_file, sep='\t', index=False, header=None)
    end = time.perf_counter()
    print(f'from {qualified_read_count} fully barcoded and with mapping quality more than {mapping_thresh_hold} reads generate {len(same_bc_reads_length_ls)} clusters, process takes {round(end - start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="Generating simulated pairs file from mapped barcode file from sprite_mapping.py")
    parser.add_argument("-in", help="mapped_barcode file", dest="input", type=str, required=True)
    parser.add_argument("-opf", help="output_file name prefix", dest="out_pf", type=str, required=True)
    parser.add_argument("-mqt", help="mapping quality thresh_hold, default value is 0", dest="mpqt", type=int, required=False)
    parser.add_argument("-lout", help="output pairs file from reads with same barcode larger than sbrt(defaul is no, to enable input yes)", dest="lp", type=str, required=False)
    parser.add_argument("-sbrlt", help="number of reads count lower threshold with the same barcode ", dest="sbrclt", type=int, required=True)
    parser.add_argument("-sbrht", help="number of reads count higher threshold with the same barcode ", dest="sbrcht", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()