#! /usr/bin/env python

import tabix
from matplotlib import pyplot as plt
import argparse
import csv
import pandas as pd
import seaborn as sns

'''
The tabix indexed input pair file should be in the format of 
['chr_1_name', 'chr_1_start', 'chr_1_pos','chr_1_strand','chr_2_name', 'chr_2_start', 'chr_2_pos','chr_2_strand','distance', 'interaction_type', 'experiment']
'interaction_type' is in the form of 'A1-A2', 'A2-A3', 'B3-B2'...
only 'A1-A1', 'A2-A2', 'B1-B1', 'B2-B2', 'B3-B3'are regarded as 'signal', others are regarded as 'noise'
'''


def run(args):
    indexed_comp_file = args.comp
    bin_file = args.bin_f
    indexed_comp_pair_file_1 = args.pair_1
    indexed_comp_pair_file_2 = args.pair_2
    cis_distance_threshold = args.dist_thr
    experiment_1_name = args.exp_1
    experiment_2_name = args.exp_2
    bins_same_comp_pair_ratio_file = args.output
    violin_figure_file = args.violin_graph

    with open(bin_file, 'r') as file1:
        with open(bins_same_comp_pair_ratio_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            pair_1_tb = tabix.open(indexed_comp_pair_file_1)
            pair_2_tb = tabix.open(indexed_comp_pair_file_2)
            compartment_tb = tabix.open(indexed_comp_file)
            bin_number = 1  # bins are counted from 1
            # the ratio column means bin_cis_signal/bin_cis_total or bin_trans_signal/bin_trans_total
            csv_writer.writerow(['bin_number', 'interaction_type', 'subcompartment_type', 'counts', 'signal_to_total_ratio', 'experiment'])
            for line in csv_reader:
                signal_type = ''
                pair_1_query_results = pair_1_tb.query(line[0], int(line[1]), int(line[2]))
                pair_2_query_results = pair_2_tb.query(line[0], int(line[1]), int(line[2]))
                comp_query_results = compartment_tb.query(line[0], int(line[1]), int(line[2]))
                try:
                    compartment = next(comp_query_results)[3]
                except:
                    pass
                if (compartment == 'A1') | (compartment == 'A2') | (compartment == 'B1') | (compartment == 'B2') | (compartment == 'B3'):
                    signal_type = compartment + '-' + compartment

                pair_1_bin_trans_signal_comp_type_dic = {'A1-A1': 0, 'A2-A2': 0, 'B1-B1': 0, 'B2-B2': 0, 'B3-B3': 0}
                pair_1_bin_cis_signal_comp_type_dic = {'A1-A1': 0, 'A2-A2': 0, 'B1-B1': 0, 'B2-B2': 0, 'B3-B3': 0}
                pair_1_bin_trans_total_count = 0
                pair_1_bin_cis_total_count = 0
                pair_2_bin_trans_signal_comp_type_dic = {'A1-A1': 0, 'A2-A2': 0, 'B1-B1': 0, 'B2-B2': 0, 'B3-B3': 0}
                pair_2_bin_cis_signal_comp_type_dic = {'A1-A1': 0, 'A2-A2': 0, 'B1-B1': 0, 'B2-B2': 0, 'B3-B3': 0}
                pair_2_bin_trans_total_count = 0
                pair_2_bin_cis_total_count = 0
                if (compartment != 'NA') & (compartment != 'B4'):
                    for result_1 in pair_1_query_results:
                        interaction_type = result_1[9]
                        if result_1[0] != result_1[4]:  # if it's a trans contact
                            if interaction_type == signal_type:  # if it's a trans signal
                                pair_1_bin_trans_signal_comp_type_dic[interaction_type] += 1
                            pair_1_bin_trans_total_count += 1
                        if result_1[0] == result_1[4]:
                            if int(float(result_1[8])) > cis_distance_threshold:  # if it's a qualified cis contact
                                if interaction_type == signal_type:  # if it's a cis signal
                                    pair_1_bin_cis_signal_comp_type_dic[interaction_type] += 1
                                pair_1_bin_cis_total_count += 1
                    for result_2 in pair_2_query_results:
                        interaction_type = result_2[9]
                        if result_2[0] != result_2[4]:  # if it's a trans contact
                            if interaction_type == signal_type:  # if it's a trans signal
                                pair_2_bin_trans_signal_comp_type_dic[interaction_type] += 1
                            pair_2_bin_trans_total_count += 1
                        if result_2[0] == result_2[4]:
                            if int(float(result_2[8])) > cis_distance_threshold:  # if it's a qualified cis contact
                                if interaction_type == signal_type:  # if it's a cis signal
                                    pair_2_bin_cis_signal_comp_type_dic[interaction_type] += 1
                                pair_2_bin_cis_total_count += 1
                    if pair_1_bin_trans_total_count != 0:
                        csv_writer.writerow([bin_number, 'trans', signal_type, pair_1_bin_trans_signal_comp_type_dic[signal_type], pair_1_bin_trans_signal_comp_type_dic[signal_type] / pair_1_bin_trans_total_count, experiment_1_name])
                    if pair_1_bin_cis_total_count != 0:
                        csv_writer.writerow([bin_number, 'cis', signal_type, pair_1_bin_cis_signal_comp_type_dic[signal_type], pair_1_bin_cis_signal_comp_type_dic[signal_type] / pair_1_bin_cis_total_count, experiment_1_name])
                    if pair_2_bin_trans_total_count != 0:
                        csv_writer.writerow([bin_number, 'trans', signal_type, pair_2_bin_trans_signal_comp_type_dic[signal_type], pair_2_bin_trans_signal_comp_type_dic[signal_type] / pair_2_bin_trans_total_count, experiment_2_name])
                    if pair_2_bin_cis_total_count != 0:
                        csv_writer.writerow([bin_number, 'cis', signal_type, pair_2_bin_cis_signal_comp_type_dic[signal_type], pair_2_bin_cis_signal_comp_type_dic[signal_type] / pair_2_bin_cis_total_count, experiment_2_name])
                bin_number += 1
    violin_plot(bins_same_comp_pair_ratio_file, violin_figure_file)


def violin_plot(input_file, output_figure_file):
    fig, ax = plt.subplots(2, 1, figsize=(6, 6), dpi=300)
    bins_same_comp_pair_ratio_df = pd.read_csv(input_file, delimiter='\t')
    cis_df = bins_same_comp_pair_ratio_df.loc[bins_same_comp_pair_ratio_df['interaction_type'] == 'cis']
    cis_df.sort_values(['subcompartment_type', 'experiment'], ascending=[True, False], inplace=True)
    trans_df = bins_same_comp_pair_ratio_df.loc[bins_same_comp_pair_ratio_df['interaction_type'] == 'trans']
    trans_df.sort_values(['subcompartment_type', 'experiment'], ascending=[True, False], inplace=True)
    v0 = sns.violinplot(x='subcompartment_type', y='signal_to_total_ratio', hue='experiment', inner="quartile", split=True, data=cis_df, ax=ax[0], palette='pastel')
    v1 = sns.violinplot(x='subcompartment_type', y='signal_to_total_ratio', hue='experiment', inner="quartile", split=True, data=trans_df, ax=ax[1], palette='pastel')
    v0.legend_.remove()
    v1.legend_.remove()
    ax[0].set_xlabel('', fontsize=15)
    ax[0].set_ylabel('cis signal to cis total ratio', fontsize=12)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    # ax[0].text(0.5, 1.2, "Cis", fontsize=12)

    ax[1].set_xlabel('signal type', fontsize=12)
    ax[1].set_ylabel('trans signal to trans total ratio', fontsize=12)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    # ax[1].text(0.5, 1.15, "Trans", fontsize=12)
    plt.setp(ax[0].collections, alpha=.6)  # this is the only way that work to change the transparency of the violin plot
    plt.setp(ax[1].collections, alpha=.6)

    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='best', frameon=False, ncol=2)
    fig.tight_layout(pad=2)
    out_put_file_name = output_figure_file + '.jpg'
    fig.savefig(out_put_file_name, format='jpg', dpi=300)


def main():
    parser = argparse.ArgumentParser(description="to calculate aggreaged breaking frequecies relative to a bunch of peak center positions")
    parser.add_argument("-ic", help="tabix indexed compartment file", dest="comp", type=str, required=True)
    parser.add_argument("-b", help="bin file without Y and M chromosomes", dest="bin_f", type=str, required=True)
    parser.add_argument("-p1", help="tabix indexed compartment tagged pair file 1", dest="pair_1", type=str, required=True)
    parser.add_argument("-p2", help="tabix indexed compartment tagged pair file 2", dest="pair_2", type=str, required=True)
    parser.add_argument("-d", help="cis distance threshold to filter out unqualified cis pairs", dest="dist_thr", type=int, required=True)
    parser.add_argument("-e1", help="experiment name for pair 1 file", dest="exp_1", type=str, required=True)
    parser.add_argument("-e2", help="experiment name for pair 2 file", dest="exp_2", type=str, required=True)
    parser.add_argument("-o", help="TF name reference file", dest="output", type=str, required=True)
    parser.add_argument("-f", help="plotted violin graph file", dest="violin_graph", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

