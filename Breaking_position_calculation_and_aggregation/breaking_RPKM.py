#! /usr/bin/env python

import argparse
import tabix
import csv
import time
import gzip


# calculate the breaking per bin per million reads
def run(args):
    start = time.perf_counter()
    bin_size = args.bins
    chrom_size_file = args.chrs
    comp_file = args.comp
    indexed_break_pos_file = args.ibpf  # breaking position file: chr_name start end
    binned_breaking_counts_file = args.output
    exp = args.exp
    tb_comp = tabix.open(comp_file)
    tb_break = tabix.open(indexed_break_pos_file)
    total_breaks = 0

    with gzip.open(indexed_break_pos_file, 'rt') as break_pos_file:
        for _ in break_pos_file:
            total_breaks += 1

    print(f'{total_breaks} found')

    chrom_dic = {}
    chrom_ls = []
    with open(chrom_size_file) as file1:
        csv_reader = csv.reader(file1, delimiter='\t')
        for line in csv_reader:
            chrom_dic.update({line[0]: int(line[1])})
            chrom_ls.append(line[0])
    print(chrom_dic)
    print(chrom_ls)

    bin_count = 1
    temp_ls = []

    with open(binned_breaking_counts_file, 'w') as file3:
        csv_writer1 = csv.writer(file3, delimiter='\t')
        for chromosome in chrom_ls:
            for i in range(0, chrom_dic[chromosome], bin_size):
                if i + bin_size < chrom_dic[chromosome]:
                    end = i + bin_size
                else:
                    end = chrom_dic[chromosome]
                try:
                    comp_result = tb_comp.query(chromosome, i, end)
                    compartment = next(comp_result)[3][:2]
                except:
                    pass
                try:
                    result = tb_break.query(chromosome, i, end)
                    for item in result:
                        temp_ls.append(item)
                    print(len(temp_ls))
                    csv_writer1.writerow([chromosome, i, end, bin_count, exp, temp_ls[0][3][:2], round(1000000*len(temp_ls)/total_breaks, 10)])
                except:
                    csv_writer1.writerow([chromosome, i, end, bin_count, exp, compartment, 0])
                temp_ls = []
                bin_count += 1
    end = time.perf_counter()
    print(f'{total_breaks} breaking position processed in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="to calculate breaking frequencies per bin per million breaks")
    parser.add_argument("-chrs", help="chromosome size file", dest="chrs", type=str, required=True)
    parser.add_argument("-bins", help="bin size in bp", dest="bins", type=int, required=True)
    parser.add_argument("-expt", help="experiment name", dest="exp", type=str, required=True)
    parser.add_argument("-ibpf", help="tabix indexed, bgzip compressed breaking position file ", dest="ibpf", type=str, required=True)
    parser.add_argument("-comp", help="tabix indexed, bgzip compressed compartment file", dest="comp", type=str, required=True)
    parser.add_argument("-out", help="aggregated list file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()