#! /usr/bin/env python

import os
import argparse
import csv
import time
import pandas as pd

'''
This code compute the cis same_compartment pair as well as the cis pair distance distance frequency,
and a
'''


def run(args):
    start = time.perf_counter()
    input_folder = args.in_file
    window = args.rolling_window
    short_dist_threshold = args.dth
    os.makedirs('./pair_statistic')  # make a folder for pie plot and data
    print('The pair statistics folder created at pair_statistic')
    os.makedirs('./cis_dist_counts')
    print('cis distance counts data folder created at cis_dist_counts')
    os.makedirs('./cis_dist_freq')
    print('The aggregated pair counts file created at cis_dist_freq')

    input_file_ls = os.listdir(input_folder)
    input_file_path_ls = []
    pair_statistic_file_path_ls = []
    same_comp_pair_statistic_file_path_ls = []
    cis_dist_counts_file_path_ls = []
    cis_dist_freq_file_path_ls = []
    same_comp_dist_freq_file_path_ls = []
    noise_dist_freq_file_path_ls = []
    cis_comp_dist_freq_file_path_ls = []

    for file in input_file_ls:
        input_file_path = input_folder + '/' + file
        input_file_path_ls.append(input_file_path)
        pair_statistic_file_path = './pair_statistic/' + file.split('.')[0] + '_pst.tsv'
        pair_statistic_file_path_ls.append(pair_statistic_file_path)
        same_comp_pair_statistic_file_path = './pair_statistic/' + file.split('.')[0] + '_scps.tsv'
        same_comp_pair_statistic_file_path_ls.append(same_comp_pair_statistic_file_path)
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
        pair_statistic_file = pair_statistic_file_path_ls[i]
        same_comp_pair_statistic_file = same_comp_pair_statistic_file_path_ls[i]
        cis_pair_file = cis_dist_counts_file_path_ls[i]
        cis_dist_freq_file = cis_dist_freq_file_path_ls[i]
        cis_comp_dist_freq_file = cis_comp_dist_freq_file_path_ls[i]
        same_comp_dist_freq_file = same_comp_dist_freq_file_path_ls[i]
        noise_dist_freq_file = noise_dist_freq_file_path_ls[i]
        experiment = get_experiment_name(input_file)

        with open(input_file, 'r') as file1:
            for _ in file1:
                total_pair_count += 1

        comp_pair_pie(input_file, pair_statistic_file, same_comp_pair_statistic_file, experiment, short_dist_threshold)
        cal_cis_dist_count(input_file, cis_pair_file)
        cal_cis_dist_frequency(cis_pair_file, total_pair_count, window, cis_dist_freq_file, cis_comp_dist_freq_file, same_comp_dist_freq_file, noise_dist_freq_file)

    new_df = pd.DataFrame(columns=['Experiment', 'Cis_Pair_Ratio', 'Trans_Pair_Ratio', 'Same_Compartment_Pair_Ratio', 'Cis_Same_Comp_to_Cis_Total', 'Trans_Same_Compartment_Pair_to_Trans_Total', f'Cis_Below_{short_dist_threshold}bp_to_Cis_Total', f'Cis_Same_Comp_Above_{short_dist_threshold}bp_to_Cis_Above_{short_dist_threshold}bp'])
    for i in range(0, len(input_file_ls)):
        pair_statistic_file = pair_statistic_file_path_ls[i]
        df = pd.read_csv(pair_statistic_file, delimiter='\t')
        new_df = pd.concat([new_df, df], sort=False)
    new_df.to_csv('./pair_statistic/total_pair_statistic.tsv', sep='\t', index=False)

    new_df1 = pd.DataFrame(columns=['Experiment', 'Interaction_Type', 'Pair_Counts', 'Percentage_to_Total_Pairs', 'Trans_pair_counts', 'Trans_Pair_to_Total_Trans_Pair_Percentage', 'Cis_Pair_Counts', 'Cis_Pair_to_Total_Cis_Pair_Percentage', f'Above_{short_dist_threshold}bp_Counts', f'Above_{short_dist_threshold}bp_Percentage', f'Below_{short_dist_threshold}bp_Counts', f'Below_{short_dist_threshold}bp_Percentage'])
    for i in range(0, len(input_file_ls)):
        same_comp_pair_statistic_file = same_comp_pair_statistic_file_path_ls[i]
        df1 = pd.read_csv(same_comp_pair_statistic_file, delimiter='\t')
        new_df1 = pd.concat([new_df1, df1], sort=False)
    new_df1.to_csv('./pair_statistic/total_same_comp_pair_statistic.tsv', sep='\t', index=False)

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


def comp_pair_pie(input_file, pair_statistic_file, same_comp_pair_statistic_file, experiment, short_dist_threshold):  # generate data for pie_plot
    comp_ls = ['A1', 'A2', 'B1', 'B2', 'B3', 'B4', 'NA']
    pair_comp_ls = []
    pair_comp_dic = {}
    cis_pair_comp_dic = {}
    pair_comp_above_threshold_dic = {}
    pair_comp_below_threshold_dic = {}
    trans_pair_comp_dic = {}
    total_comp_pair_counts = 0
    same_comp_pair_counts = 0
    trans_pair_counts = 0
    same_comp_trans_pair_counts = 0
    same_comp_cis_pair_counts = 0
    above_threshold_same_cis_comp_pair_count = 0
    cis_pair_counts = 0
    above_threshold_cis_pair_count = 0
    below_threshold_cis_pair_count = 0
    for i in range(0, len(comp_ls)):   # all possible interaction type: A1-A1, A1-A2...A2-A2..A2-B1,A2-B2...B1-B1
        for j in range(i, len(comp_ls)):
            pair_comp_ls.append(comp_ls[i] + '-' + comp_ls[j])
            pair_comp_dic.update({comp_ls[i] + '-' + comp_ls[j]: 0})  # make a dictionary counter
            cis_pair_comp_dic.update({comp_ls[i] + '-' + comp_ls[j]: 0})
            pair_comp_above_threshold_dic.update({comp_ls[i] + '-' + comp_ls[j]: 0})
            pair_comp_below_threshold_dic.update({comp_ls[i] + '-' + comp_ls[j]: 0})
            trans_pair_comp_dic.update({comp_ls[i] + '-' + comp_ls[j]: 0})
    with open(input_file, 'r') as file1:
        with open(pair_statistic_file, 'w') as file2:
            with open(same_comp_pair_statistic_file, 'w') as file3:
                csv_reader = csv.reader(file1, delimiter='\t')
                csv_writer = csv.writer(file2, delimiter='\t')
                csv_writer1 = csv.writer(file3, delimiter='\t')
                for pair in csv_reader:
                    try:
                        pair_comp_dic[pair[7]] += 1  # This line count the total pairs that belong to certain comp_pair type, pair[7] is the interaction type
                    except KeyError:
                        reversed_pair_comp_type = pair[7].split('-')[1] + '-' + pair[7].split('-')[0]
                        pair_comp_dic[reversed_pair_comp_type] += 1  # This considers cases such A1-B1, B1-A1, both hit 'A1-B1' key

                    if pair[0] == pair[3]:  # if it is a cis pair
                        cis_pair_counts += 1
                        try:
                            cis_pair_comp_dic[pair[7]] += 1
                        except KeyError:
                            reversed_pair_comp_type = pair[7].split('-')[1] + '-' + pair[7].split('-')[0]
                            cis_pair_comp_dic[reversed_pair_comp_type] += 1
                        if int(pair[6]) > short_dist_threshold:
                            above_threshold_cis_pair_count += 1
                            try:
                                pair_comp_above_threshold_dic[pair[7]] += 1
                            except KeyError:
                                reversed_pair_comp_type = pair[7].split('-')[1] + '-' + pair[7].split('-')[0]
                                pair_comp_above_threshold_dic[reversed_pair_comp_type] += 1
                        elif int(pair[6]) <= short_dist_threshold:
                            below_threshold_cis_pair_count += 1
                            try:
                                pair_comp_below_threshold_dic[pair[7]] += 1
                            except KeyError:
                                reversed_pair_comp_type = pair[7].split('-')[1] + '-' + pair[7].split('-')[0]
                                pair_comp_below_threshold_dic[reversed_pair_comp_type] += 1

                    if pair[0] != pair[3]:  # if it is not a cis pair
                        trans_pair_counts += 1
                        try:
                            trans_pair_comp_dic[pair[7]] += 1
                        except KeyError:
                            reversed_pair_comp_type = pair[7].split('-')[1] + '-' + pair[7].split('-')[0]
                            trans_pair_comp_dic[reversed_pair_comp_type] += 1

                    if pair[7] in ['A1-A1', 'A2-A2', 'B1-B1', 'B2-B2', 'B3-B3', 'B4-B4']:  # if it is a same_comp pair
                        same_comp_pair_counts += 1
                        if pair[0] != pair[3]:  # if it's a trans_same_comp pair
                            same_comp_trans_pair_counts += 1
                        elif pair[0] == pair[3]:  # if it's a cis_same_comp pair
                            same_comp_cis_pair_counts += 1
                            if int(pair[6]) > short_dist_threshold:
                                above_threshold_same_cis_comp_pair_count += 1

                    total_comp_pair_counts += 1  # This line count the total number of pairs parsed

                csv_writer.writerow(['Experiment', 'Cis_Pair_Ratio', 'Trans_Pair_Ratio', 'Same_Compartment_Pair_Ratio', 'Cis_Same_Comp_to_Cis_Total', 'Trans_Same_Compartment_Pair_to_Trans_Total', f'Cis_Below_{short_dist_threshold}bp_to_Cis_Total', f'Cis_Same_Comp_Above_{short_dist_threshold}bp_to_Cis_Above_{short_dist_threshold}bp'])  # write the header row
                csv_writer.writerow([experiment, round(100*cis_pair_counts/total_comp_pair_counts, 5), round(100*trans_pair_counts/total_comp_pair_counts, 5), round(100*same_comp_pair_counts/total_comp_pair_counts, 5), round(100*same_comp_cis_pair_counts/cis_pair_counts, 5), round(100*same_comp_trans_pair_counts/trans_pair_counts, 5), round(100*below_threshold_cis_pair_count/cis_pair_counts, 5), round(100*above_threshold_same_cis_comp_pair_count/above_threshold_cis_pair_count, 5)])
                csv_writer1.writerow(['Experiment', 'Interaction_Type', 'Pair_Counts', 'Percentage_to_Total_Pairs', 'Trans_pair_counts', 'Trans_Pair_to_Total_Trans_Pair_Percentage', 'Cis_Pair_Counts', 'Cis_Pair_to_Total_Cis_Pair_Percentage', f'Above_{short_dist_threshold}bp_Counts', f'Above_{short_dist_threshold}bp_Percentage', f'Below_{short_dist_threshold}bp_Counts', f'Below_{short_dist_threshold}bp_Percentage'])
                for pair_comp in pair_comp_ls:
                    csv_writer1.writerow([experiment, pair_comp, pair_comp_dic[pair_comp], round(100*pair_comp_dic[pair_comp]/total_comp_pair_counts, 5), trans_pair_comp_dic[pair_comp], round(100*trans_pair_comp_dic[pair_comp]/trans_pair_counts, 5), cis_pair_comp_dic[pair_comp], round(100*cis_pair_comp_dic[pair_comp]/cis_pair_counts, 5), pair_comp_above_threshold_dic[pair_comp], round(100*pair_comp_above_threshold_dic[pair_comp]/above_threshold_cis_pair_count, 5), pair_comp_below_threshold_dic[pair_comp], round(100*pair_comp_below_threshold_dic[pair_comp]/below_threshold_cis_pair_count, 5)])


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
    parser.add_argument("-dt", help="the distance threshold (in bp) above which the pair being separated for same compartment pair analysis", dest="dth",
                        type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
