#! /usr/bin/env python

import cooler
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import cm


def run(args):
    cool_file_1 = args.cool_1
    cool_file_2 = args.cool_2
    bal_choice = args.blc
    fig_file = args.fig

    c1 = cooler.Cooler(cool_file_1)
    c2 = cooler.Cooler(cool_file_2)
    bins = c1.bins()[:]
    bin_dic = {}
    chrom_bin_dic = {}
    for index, row in bins.iterrows():  # bins is a pandas dataframe
        bin_dic[index] = [row[0], row[1], row[2]]
    chromosome_ls = []
    for i in range(1, 23):  # generate a list that contains all the chromosomes' names that the compartment file has
        chromosome_ls.append('chr' + str(i))
    chromosome_ls.append('chrX')
    # chromosome_ls.append('chrY')
    # chromosome_ls.append('chrM')
    for chromosome in chromosome_ls:
        chrom_bins = bins.loc[bins['chrom'] == chromosome]
        chrom_bin_dic[chromosome] = (chrom_bins.index[0], chrom_bins.index[-1])

    chr_label_ls = []
    chr_border_coord_ls = []
    for chromosome in chromosome_ls:
        label_pos = int((chrom_bin_dic[chromosome][1] + chrom_bin_dic[chromosome][0]) / 2)
        chr_label_ls.append(label_pos)
        chr_border_coord_ls.append(chrom_bin_dic[chromosome][0])
        chr_border_coord_ls.append(chrom_bin_dic[chromosome][1])

    if bal_choice == "y":
        bal_choice = True
    elif bal_choice == "n":
        bal_choice = False
    else:
        print("please enter y or n for balancing choice")
    c1_matrix = c1.matrix(balance=bal_choice, sparse=False)[:]
    c2_matrix = c2.matrix(balance=bal_choice, sparse=False)[:]
    div_matrix = np.divide(c1_matrix, c2_matrix)
    div_matrix_noinf = np.where(np.isinf(div_matrix), np.nan, div_matrix)
    new_2d_ls = []
    new_ls = []
    for i in range(0, chrom_bin_dic['chrX'][
        1]):  # chrX is the third to the last chromosome in the bin list, chrom_bin_dic['chrX'][1] is the chrX bin end number
        for j in range(0, chrom_bin_dic['chrX'][1]):
            new_ls.append(div_matrix_noinf[i, j])
            if len(new_ls) >= chrom_bin_dic['chrX'][1]:
                new_2d_ls.append(new_ls)
                new_ls = []
    div_matrix_noinf_trimmed = np.array(new_2d_ls)
    div_matrix_noinf_trimmed_log_transformed = np.log10(div_matrix_noinf_trimmed)
    div_matrix_noinf_trimmed_log_transformed_noinf = np.where(np.isinf(div_matrix_noinf_trimmed_log_transformed), np.nan, div_matrix_noinf_trimmed_log_transformed)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(16, 8), dpi=300)
    fig.tight_layout(pad=5)  # make enough space between subplots
    color_map = plt.cm.get_cmap('bwr')  # blue white and red from value low to high
    max_value_1 = np.nanmax(div_matrix_noinf_trimmed)
    min_value_1 = np.nanmin(div_matrix_noinf_trimmed)
    print(f"highest value of the division matrix is {max_value_1}")
    print(f"lowest value of the division matrix is {min_value_1}")

    norm_1 = cm.colors.DivergingNorm(vmax=max_value_1, vcenter=1, vmin=min_value_1)  # 1 is white
    im_1 = ax[0].matshow(div_matrix_noinf_trimmed, cmap=color_map, norm=norm_1)
    cax_1 = fig.add_axes([0.47, 0.71, 0.01, 0.2])  # [left, bottom, width, height] of the new ax to plot the color bar
    cbar_1 = plt.colorbar(im_1, cax=cax_1)  # cax is the ax parameter defining where colorbar will be drawn.
    cbar_1.ax.tick_params(labelsize=5)  # change the color bar text size
    ax[0].set_title("without log transformation", fontsize=10, y=1, pad=8)
    ax[0].set_xticks(chr_label_ls)
    ax[0].set_xticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'x'], size=6)
    ax[0].set_xticks(chr_border_coord_ls, minor=True)  # set_xticks is not required for correct plotting. I think there's a bug
    ax[0].tick_params(which='minor', length=1, width=0.2, bottom=False, right=False)

    ax[0].set_yticks(chr_label_ls)
    ax[0].set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'x'], size=6)
    ax[0].set_yticks(chr_border_coord_ls, minor=True)  # ax.set_yticks is required for correct plotting, I think there's a bug
    ax[0].tick_params(which='minor', length=1, width=0.2, bottom=False, right=False)  # set minor tick size and location
    ax[0].tick_params(axis='both', which='major', bottom=False, top=False, labelbottom=False, right=False, left=False)  # inactivate major ticks
    # ax[0].tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False)
    ax[0].spines['right'].set_visible(False)  # these are used to remove the black frame around the plot.
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)

    max_value_2 = np.nanmax(div_matrix_noinf_trimmed_log_transformed_noinf)
    min_value_2 = np.nanmin(div_matrix_noinf_trimmed_log_transformed_noinf)
    print(f"highest value of log10 transformed division matrix is {max_value_2}")
    print(f"lowest value of log10 transformed division matrix is {min_value_2}")
    norm_2 = cm.colors.DivergingNorm(vmax=max_value_2, vcenter=0)
    im_2 = ax[1].matshow(div_matrix_noinf_trimmed_log_transformed_noinf, cmap=color_map, norm=norm_2, label="with_log_transformation")
    cax_2 = fig.add_axes([0.95, 0.71, 0.01, 0.2])  # [left, bottom, width, height] of the new axes
    cbar_2 = plt.colorbar(im_2, cax=cax_2)  # cax is the ax parameter defining where colorbar will be drawn.
    cbar_2.ax.tick_params(labelsize=5)  # change the color bar text size
    ax[1].set_title("with log transformation", fontsize=10, y=1, pad=8)
    ax[1].set_xticks(chr_label_ls)
    ax[1].set_xticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X'], size=6)
    ax[1].set_xticks(chr_border_coord_ls, minor=True)
    ax[1].set_yticks(chr_label_ls)
    ax[1].set_yticklabels([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 'X'], size=6)
    ax[1].set_yticks(chr_border_coord_ls, minor=True)
    ax[1].tick_params(which='minor', length=1, width=0.2, bottom=False, right=False)  # set minor tick size and location
    ax[1].tick_params(axis='both', which='major', bottom=False, top=False, labelbottom=False, right=False, left=False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['left'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)

    fig.savefig(fig_file, format='png', dpi=300)


def main():
    parser = argparse.ArgumentParser(description="draw division matrix from cool file with blue white red color bar")
    parser.add_argument("-c1", help="cool file on top", dest="cool_1", type=str, required=True)
    parser.add_argument("-c2", help="cool file on bottom", dest="cool_2", type=str, required=True)
    parser.add_argument("-bc", help="balance or not, enter y or n", dest="blc", type=str, required=True)
    parser.add_argument("-f", help="output figure file name", dest="fig", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
