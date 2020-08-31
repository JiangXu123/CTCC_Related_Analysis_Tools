#! /usr/bin/env python

import csv
import time
import argparse
import ast
import pandas as pd


def run(args):
    same_bc_read_list_file = args.input
    output_file_prefix = args.out_pf
    same_bc_reads_count_lower_thresh_hold = args.sbrclt
    same_bc_reads_count_higher_thresh_hold = args.sbrcht
    smaller_output_file = output_file_prefix + "sbrlt_" + str(same_bc_reads_count_lower_thresh_hold) + "_l.pair"
    larger_output_file = output_file_prefix + "sbrlt_" + str(same_bc_reads_count_lower_thresh_hold) + "_sbrht_" + str(same_bc_reads_count_higher_thresh_hold) + ".pair"
    same_bc_read_length_list_file = output_file_prefix + "_sbcrl.tsv"
    mapping_thresh_hold = args.mpqt
    output_distance_larger_than_same_bc_read_threshold = args.lp
    item_count = 0
    same_bc_reads_length_ls = []

    start = time.perf_counter()
    with open(same_bc_read_list_file, 'r') as csv_file3:
        csv_reader3 = csv.reader(csv_file3, delimiter='\t')
        with open(smaller_output_file, 'w') as csv_file4:
            csv_writer2 = csv.writer(csv_file4, delimiter='\t')
            with open(larger_output_file,'w') as csv_file5:
                csv_writer3 = csv.writer(csv_file5, delimiter='\t')
                for item in csv_reader3:
                    temp_ls = []
                    for indiv_ls in item:
                        if ast.literal_eval(indiv_ls) not in temp_ls:
                            temp_ls.append(ast.literal_eval(indiv_ls))  # convert each item in the list from a string to a list data format
                same_bc_read_length = len(temp_ls)
                same_bc_reads_length_ls.append(same_bc_read_length)
                if same_bc_read_length <= same_bc_reads_count_lower_thresh_hold:
                    for a in range(0, len(temp_ls) - 1):
                        for b in range(a + 1, len(smaller_temp_ls)):
                            if temp_ls[a][0] == temp_ls[b][0]:
                                csv_writer2.writerow([temp_ls[a][0], temp_ls[a][1] - temp_ls[b][1], 1])
                if output_distance_larger_than_same_bc_read_threshold == 'yes':
                    if (same_bc_read_length > same_bc_reads_count_lower_thresh_hold) & (same_bc_read_length <= same_bc_reads_count_higher_thresh_hold):
                        for a in range(0, len(temp_ls) - 1):
                            for b in range(a + 1, len(smaller_temp_ls)):
                                if temp_ls[a][0] == temp_ls[b][0]:
                                    csv_writer3.writerow([temp_ls[a][0], temp_ls[a][1] - temp_ls[b][1], 1])
                item_count += 1
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