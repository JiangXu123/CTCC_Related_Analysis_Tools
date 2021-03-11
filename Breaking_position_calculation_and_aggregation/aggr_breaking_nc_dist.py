#! /usr/bin/env python

import argparse
import csv
import time


# Aggregate breaking position relative to peak position
def run(args):
    start = time.perf_counter()
    break_distance_file = args.bdf
    break_table_file = args.table
    expt = args.experiment
    peak_sig_height = args.psh
    A1_dic = {}
    A2_dic = {}
    B1_dic = {}
    B2_dic = {}
    B3_dic = {}
    B4_dic = {}
    A1_bin_ls = []
    A2_bin_ls = []
    B1_bin_ls = []
    B2_bin_ls = []
    B3_bin_ls = []
    B4_bin_ls = []

    compartment_ls = ['A1', 'A2', 'B1', 'B2', 'B3', 'B4']
    line_number = 0

    for i in range(-600, 601):
        A1_dic.update({i: 0})
        A2_dic.update({i: 0})
        B1_dic.update({i: 0})
        B2_dic.update({i: 0})
        B3_dic.update({i: 0})
        B4_dic.update({i: 0})

    with open(break_distance_file, 'r') as break_pos:
        with open(break_table_file, 'w') as aggregated_data:
            csv_reader = csv.reader(break_pos, delimiter='\t')
            csv_writer = csv.writer(aggregated_data, delimiter='\t')
            for line in csv_reader:
                if line[2] == 'A1':
                    A1_dic[int(line[4])] += 1  # line[4] records the distance from nucleosome center.
                    A1_bin_ls.append(line[1])  # line[1] records the nucleosome ID number(each nucleosome has an ID number)
                elif line[2] == 'A2':
                    A2_dic[int(line[4])] += 1
                    A2_bin_ls.append(line[1])
                elif line[2] == 'B1':
                    B1_dic[int(line[4])] += 1
                    B1_bin_ls.append(line[1])
                elif line[2] == 'B2':
                    B2_dic[int(line[4])] += 1
                    B2_bin_ls.append(line[1])
                elif line[2] == 'B3':
                    B3_dic[int(line[4])] += 1
                    B3_bin_ls.append(line[1])
                elif line[2] == 'B4':
                    B4_dic[int(line[4])] += 1
                    B4_bin_ls.append(line[1])
                line_number += 1
                if line_number%1000 == 0:
                    print(f"{line_number} of lines processed")
            A1_bin_count = len(list(dict.fromkeys(A1_bin_ls)))
            A2_bin_count = len(list(dict.fromkeys(A2_bin_ls)))
            B1_bin_count = len(list(dict.fromkeys(B1_bin_ls)))
            B2_bin_count = len(list(dict.fromkeys(B2_bin_ls)))
            B3_bin_count = len(list(dict.fromkeys(B3_bin_ls)))
            B4_bin_count = len(list(dict.fromkeys(B4_bin_ls)))
            print(f'bin belong to A1 have {A1_bin_count}')
            print(f'bin belong to A2 have {A2_bin_count}')
            print(f'bin belong to B1 have {B1_bin_count}')
            print(f'bin belong to B2 have {B2_bin_count}')
            print(f'bin belong to B3 have {B3_bin_count}')
            print(f'bin belong to B4 have {B4_bin_count}')
            for compartment in compartment_ls:
                for i in range(-600, 601):    # i is the relative distance to the center of the peak in bp
                    if compartment == 'A1':
                        csv_writer.writerow([expt, compartment, i, A1_dic[i], A1_bin_count, line_number, 1000000*A1_dic[i]/(A1_bin_count * line_number), peak_sig_height])
                    if compartment == 'A2':
                        csv_writer.writerow([expt, compartment, i, A2_dic[i], A2_bin_count, line_number, 1000000*A2_dic[i]/(A2_bin_count * line_number), peak_sig_height])
                    if compartment == 'B1':
                        csv_writer.writerow([expt, compartment, i, B1_dic[i], B1_bin_count, line_number, 1000000*B1_dic[i]/(B1_bin_count * line_number), peak_sig_height])
                    if compartment == 'B2':
                        csv_writer.writerow([expt, compartment, i, B2_dic[i], B2_bin_count, line_number, 1000000*B2_dic[i]/(B2_bin_count * line_number), peak_sig_height])
                    if compartment == 'B3':
                        csv_writer.writerow([expt, compartment, i, B3_dic[i], B3_bin_count, line_number, 1000000*B3_dic[i]/(B3_bin_count * line_number), peak_sig_height])
                    if compartment == 'B4':
                        csv_writer.writerow([expt, compartment, i, B4_dic[i], B4_bin_count, line_number, 1000000*B4_dic[i]/(B4_bin_count * line_number), peak_sig_height])
            # column names are: experiment_type, compartmet, position_to_nc_center, break_count, compartment_bin_count,
            # total_break_count, break_count/bin_count*total_break_count

    end = time.perf_counter()
    print(f'{line_number} breaking position processed in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="to calculate breaking frequecy relative to a buch of peaks center position")
    parser.add_argument("-bdf", help="breaking distance to peak center file ", dest="bdf", type=str, required=True)
    parser.add_argument("-out", help="aggregated output table file ", dest="table", type=str, required=True)
    parser.add_argument("-expt", help="experiment name", dest="experiment", type=str, required=True)
    parser.add_argument("-psh", help="Peak signal height", dest="psh", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()