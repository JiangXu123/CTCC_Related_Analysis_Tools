#! /usr/bin/env python

import os
import argparse
import csv
import time
import pandas as pd


def run(args):
    start = time.perf_counter()
    input_folder = args.in_file
    window = args.rolling_window
    os.makedirs('./pie_plot')  # make a folder for pie plot and data
    print('The pie plot folder created at pie_plot')
    os.makedirs('./cis_dist_counts')
    print('cis distance counts data folder created at cis_dist_counts')
    os.makedirs('./cis_dist_freq')
    print('The aggregated pair counts file created at cis_dist_freq')

    input_file_ls = os.listdir(input_folder)
    input_file_path_ls = []
    pie_plot_file_path_ls = []
    cis_dist_counts_file_path_ls = []
    cis_dist_freq_file_path_ls = []
    same_comp_dist_freq_file_path_ls = []
    noise_dist_freq_file_path_ls = []
    cis_comp_dist_freq_file_path_ls = []

    for file in input_file_ls:
        input_file_path = input_folder + '/' + file
        input_file_path_ls.append(input_file_path)
        pie_plot_file_path = './pie_plot/' + file.split('.')[0] + '_pie.tsv'
        pie_plot_file_path_ls.append(pie_plot_file_path)
        cis_dist_counts_file_path = './cis_dist_counts/' + file.split('.')[0] + '_cis_dc.tsv'
        cis_dist_counts_file_path_ls.append(cis_dist_counts_file_path)
        cis_dist_freq_file_path = './cis_dist_freq/' + file.split('.')[0] + '_cis_df.tsv'
        cis_dist_freq_file_path_ls.append(cis_dist_freq_file_path)
        cis_comp_dist_freq_file_path = './cis_dist_freq/' + file.split('.')[0] + '_cis_comp_df.tsv'
        cis_comp_dist_freq_file_path_ls.append(cis_comp_dist_freq_file_path)
        same_comp_dist_freq_file_path = './cis_dist_freq/' + file.split('.')[0] + '_sc_df.tsv'
        same_comp_dist_freq_file_path_ls.append(same_comp_dist_freq_file_path)
        noise_dist_freq_file_path = './cis_dist_freq/' + file.split('.')[0] + '_noise_df.tsv'
        noise_dist_freq_file_path_ls.append(noise_dist_freq_file_path)

    for i in range(0, len(input_file_ls)):
        total_pair_count = 0
        input_file = input_file_path_ls[i]
        pie_plot_file = pie_plot_file_path_ls[i]
        cis_pair_file = cis_dist_counts_file_path_ls[i]
        cis_dist_freq_file = cis_dist_freq_file_path_ls[i]
        cis_comp_dist_freq_file = cis_comp_dist_freq_file_path_ls[i]
        same_comp_dist_freq_file = same_comp_dist_freq_file_path_ls[i]
        noise_dist_freq_file = noise_dist_freq_file_path_ls[i]
        experiment = get_experiment_name(input_file)

        with open(input_file, 'r') as file1:
            for _ in file1:
                total_pair_count += 1

        comp_pair_pie(input_file, pie_plot_file, experiment)
        cal_cis_dist_count(input_file, cis_pair_file)
        cal_cis_dist_frequency(cis_pair_file, total_pair_count, window, cis_dist_freq_file, cis_comp_dist_freq_file, same_comp_dist_freq_file, noise_dist_freq_file)

    new_df1 = pd.DataFrame(columns=['Pair_Compartment', 'Counts', 'Counts_Per_Million', 'Same_Compartment_Pair_Fraction', 'Experiment'])
    for i in range(0, len(input_file_ls)):
        pie_plot_file = pie_plot_file_path_ls[i]
        df1 = pd.read_csv(pie_plot_file, delimiter='\t')
        new_df1 = pd.concat([new_df1, df1], sort=False)
    new_df1.to_csv('./pie_plot/total_pair_comp_statistic.tsv', sep='\t', index=False)

    new_df2 = pd.DataFrame(columns=['Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million'])
    for i in range(0, len(input_file_ls)):
        cis_dist_freq_file = cis_dist_freq_file_path_ls[i]
        df2 = pd.read_csv(cis_dist_freq_file, delimiter='\t')
        new_df2 = pd.concat([new_df2, df2], sort=False)
    new_df2.to_csv('./cis_dist_freq/total_cis_dist_freq.tsv', sep='\t', index=False)

    new_df3 = pd.DataFrame(columns=['Compartment_Type', 'Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million'])
    for i in range(0, len(input_file_ls)):
        cis_comp_dist_freq_file = cis_comp_dist_freq_file_path_ls[i]
        df3 = pd.read_csv(cis_comp_dist_freq_file, delimiter='\t')
        new_df3 = pd.concat([new_df3, df3], sort=False)
    new_df3.to_csv('./cis_dist_freq/total_comp_dist_freq.tsv', sep='\t', index=False)

    new_df4 = pd.DataFrame(columns=['Compartment_Type', 'Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million'])
    for i in range(0, len(input_file_ls)):
        same_comp_dist_freq_file = same_comp_dist_freq_file_path_ls[i]
        df4 = pd.read_csv(same_comp_dist_freq_file, delimiter='\t')
        new_df4 = pd.concat([new_df4, df4], sort=False)
    new_df4.to_csv('./cis_dist_freq/total_same_comp_dist_freq.tsv', sep='\t', index=False)

    new_df5 = pd.DataFrame(columns=['Compartment_Type', 'Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million'])
    for i in range(0, len(input_file_ls)):
        noise_dist_freq_file = noise_dist_freq_file_path_ls[i]
        df5 = pd.read_csv(noise_dist_freq_file, delimiter='\t')
        new_df5 = pd.concat([new_df5, df5], sort=False)
    new_df5.to_csv('./cis_dist_freq/total_noise_dist_freq.tsv', sep='\t', index=False)

    end = time.perf_counter()
    print(f'files processed in {round(end - start, 2)} seconds')


def get_experiment_name(input_file):
    with open(input_file, 'r') as file1:
        csv_reader = csv.reader(file1, delimiter='\t')
        experiment = next(csv_reader)[8]
        return experiment


def comp_pair_pie(input_file, same_comp_pie_file, experiment):  # generate data for pie_plot
    comp_ls = ['A1', 'A2', 'B1', 'B2', 'B3', 'B4', 'NA']
    pair_comp_ls = []
    pair_comp_dic = {}
    total_comp_pair_counts = 0
    same_comp_pair_counts = 0
    for i in range(0, len(comp_ls)):
        for j in range(i, len(comp_ls)):
            pair_comp_ls.append(comp_ls[i] + '-' + comp_ls[j])
            pair_comp_dic.update({comp_ls[i] + '-' + comp_ls[j]: 0})

    with open(input_file, 'r') as file1:
        with open(same_comp_pie_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            for pair in csv_reader:
                try:
                    pair_comp_dic[pair[7]] += 1  # This line count the total pairs that belong to certain comp_pair type, pair[7] is the interaction type
                except KeyError:
                    reversed_pair_comp_type = pair[7].split('-')[1] + '-' + pair[7].split('-')[0]
                    pair_comp_dic[reversed_pair_comp_type] += 1  # This will consider cases like A1-B1, B1-A1, both in the 'A1-B1' key
                if pair[7] in ['A1-A1', 'A2-A2', 'B1-B1', 'B2-B2', 'B3-B3', 'B4-B4']:
                    same_comp_pair_counts += 1
                total_comp_pair_counts += 1  # This line count the total number compartmented pairs
            csv_writer.writerow(['Pair_Compartment', 'Counts', 'Counts_Per_Million', 'Same_Compartment_Pair_Fraction', 'Experiment'])
            for pair_comp in pair_comp_ls:
                csv_writer.writerow([pair_comp, pair_comp_dic[pair_comp], 1000000*pair_comp_dic[pair_comp]/total_comp_pair_counts, 100*same_comp_pair_counts/total_comp_pair_counts, experiment])
         # columns are ['Pair_Compartment', 'Counts', 'Counts_Per_Million', 'Same_Compartment_Pair_Fraction', 'Experiment']


def cal_cis_dist_count(input_file, cis_pair_file):
    with open(input_file, 'r') as file1:
        with open(cis_pair_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            csv_writer.writerow(['Chromosome', 'Distance', 'Compartment_Type', 'Experiment'])
            for pair in csv_reader:
                cis_pair_count = 0
                if pair[0] == pair[3]:  # if it is a cis contact
                    csv_writer.writerow([pair[0], pair[6], pair[7], pair[8]])  # distance, comp_type, experiment
                    cis_pair_count += 1


def cal_cis_dist_frequency(cis_pair_file, total_pair_count, window, cis_dist_freq_file, cis_comp_dist_freq_file, same_comp_dist_freq_file, noise_dist_freq_file):
    df1 = pd.read_csv(cis_pair_file, delimiter='\t')
    df1['Counts'] = 1

    df1_freq = df1.groupby(['Distance', 'Experiment']).count()['Counts'].reset_index()
    df1_freq['Counts_Per_Million'] = 1000000 * df1_freq['Counts'] / total_pair_count
    df1_freq['Average_Counts_Per_Million'] = df1_freq['Counts_Per_Million'].rolling(window).mean().reset_index(level=0, drop=True)
    df1_freq[['Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million']].to_csv(cis_dist_freq_file, sep='\t', index=False)

    df1_comp_freq = df1.groupby(['Compartment_Type', 'Distance', 'Experiment']).count()['Counts'].reset_index()
    df1_comp_freq['Counts_Per_Million'] = 1000000 * df1_comp_freq['Counts'] / total_pair_count
    df1_comp_freq['Average_Counts_Per_Million'] = df1_comp_freq['Counts_Per_Million'].rolling(window).mean().reset_index(level=0, drop=True)
    df1_comp_freq[['Compartment_Type', 'Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million']].to_csv(cis_comp_dist_freq_file, sep='\t', index=False)

    df1_same_comp_freq = df1.loc[(df1['Compartment_Type'] == 'A1-A1') | (df1['Compartment_Type'] == 'A2-A2') | (df1['Compartment_Type'] == 'B1-B1') | (df1['Compartment_Type'] == 'B2-B2') | (df1['Compartment_Type'] == 'B3-B3') | (df1['Compartment_Type'] == 'B4-B4')]
    df1_same_comp_freq = df1_same_comp_freq.groupby(['Compartment_Type', 'Distance', 'Experiment']).count()['Counts'].reset_index()
    df1_same_comp_freq['Counts_Per_Million'] = 1000000 * df1_same_comp_freq['Counts'] / total_pair_count
    df1_same_comp_freq['Average_Counts_Per_Million'] = df1_same_comp_freq['Counts_Per_Million'].rolling(window).mean().reset_index(level=0, drop=True)
    df1_same_comp_freq[['Compartment_Type', 'Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million']].to_csv(same_comp_dist_freq_file, sep='\t', index=False)

    df1_noise_dist_freq = df1.loc[(df1['Compartment_Type'] != 'A1-A1') & (df1['Compartment_Type'] != 'A2-A2') & (df1['Compartment_Type'] != 'B1-B1') & (df1['Compartment_Type'] != 'B2-B2') & (df1['Compartment_Type'] != 'B3-B3') & (df1['Compartment_Type'] != 'B4-B4')]
    df1_noise_dist_freq = df1_noise_dist_freq.groupby(['Compartment_Type', 'Distance', 'Experiment']).count()['Counts'].reset_index()
    df1_noise_dist_freq['Counts_Per_Million'] = 1000000 * df1_noise_dist_freq['Counts'] / total_pair_count
    df1_noise_dist_freq['Average_Counts_Per_Million'] = df1_noise_dist_freq['Counts_Per_Million'].rolling(window).mean().reset_index(level=0, drop=True)
    df1_noise_dist_freq[['Compartment_Type', 'Distance', 'Experiment', 'Counts', 'Counts_Per_Million', 'Average_Counts_Per_Million']].to_csv(noise_dist_freq_file, sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description="aggregated and calculate the HiC pair distance frequency")
    parser.add_argument("-i", help="input folder", dest="in_file", type=str, required=True)
    parser.add_argument("-w", help="rolling average window size for plotting dist-freq curve", dest="rolling_window", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
