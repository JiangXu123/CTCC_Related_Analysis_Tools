#! /usr/bin/env python

import argparse
import csv
import time
import pandas as pd

'''
This code is used to compute and plot cis pair distance-pair-frequency for plotting
'''


def run(args):
    start = time.perf_counter()
    input_pair_file = args.pair
    window = args.window_s
    exp_name = args.exp_n
    cis_dist_file = args.cis_dist
    cis_dist_freq_file = args.cis_dist_freq
    cal_cis_dist_frequency(input_pair_file, exp_name, window, cis_dist_file, cis_dist_freq_file)
    end = time.perf_counter()
    print(f'files processed in {round(end - start, 2)} seconds')


def cal_cis_dist_frequency(pair_file, exp_name, window, cis_dist_file, cis_dist_freq_file):
    with open(pair_file, 'r') as file1:
        with open(cis_dist_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            csv_writer.writerow(['experiment', 'chr_name', 'distance', 'count'])
            total_pair_count = 0
            for line in csv_reader:
                if line[1] == line[4]: # if it's a cis pair
                    if int(line[2]) >= int(line[5]):
                        cis_distance = int(line[2]) - int(line[5])
                    elif int(line[2]) < int(line[5]):
                        cis_distance = int(line[5]) - int(line[2])
                    csv_writer.writerow([exp_name, line[1], cis_distance, 1])  # the 'Count' column writen to 1, should've be unnecessary
                    total_pair_count += 1
    cis_dist_df = pd.read_csv(cis_dist_file, delimiter='\t')
    cis_dist_freq_df = cis_dist_df.groupby(['distance', 'experiment']).count()['count'].reset_index()
    cis_dist_freq_df['Counts_Per_Million'] = 1000000 * cis_dist_freq_df['count'] / total_pair_count
    cis_dist_freq_df['Average_Counts_Per_Million'] = cis_dist_freq_df['Counts_Per_Million'].rolling(window).mean().reset_index(level=0, drop=True)
    cis_dist_freq_df[['distance', 'experiment', 'count', 'Counts_Per_Million', 'Average_Counts_Per_Million']].to_csv(cis_dist_freq_file, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description="aggregated and calculate the HiC pair distance frequency")
    parser.add_argument("-i", help="input pair file", dest="pair", type=str, required=True)
    parser.add_argument("-w", help="rolling average window size for plotting dist-freq curve", dest="window_s", type=int, required=True)
    parser.add_argument("-n", help="experiment name", dest="exp_n", type=str, required=True)
    parser.add_argument("-d", help="output cis distance", dest="cis_dist", type=str, required=True)
    parser.add_argument("-df", help="output cis distance frequency file", dest="cis_dist_freq", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
