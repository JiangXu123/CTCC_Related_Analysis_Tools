#! /usr/bin/env python

import csv
import time
from operator import itemgetter
import argparse


def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    data_ls = []
    with open(input_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for line in csv_reader:
            if (int(line[4]) > 101010101) & (line[4][-2:] != '00') & (line[4][-4:-2] != '00') & (
                    line[4][-6:-4] != '00') & (line[4][-8:-6] != '00') & (line[4][0:-8] != '00') & (line[1] != ''):
                data_ls.append(line)
    sorted_ls = sorted(data_ls, key=itemgetter(4))

    ls = []
    all_ls = []
    temp = sorted_ls[0]

    ls.append(temp)
    for i in range(1, len(sorted_ls)):
        if temp[4] == sorted_ls[i][4]:
            ls.append(sorted_ls[i])
        else:
            if len(ls) > 1:
                all_ls.append(ls)
            temp = sorted_ls[i]
            ls = []
            ls.append(temp)

    pair_ls = []
    for item in all_ls:
        for a in range(0, len(item) - 1):
            for b in range(a + 1, len(item)):
                pair_ls.append([item[a][1], item[a][2], item[b][1], item[b][2], 1])
    pair_ls.sort(key=lambda x: [x[0], x[1], x[2], x[3]])
    with open(output_file, 'w') as new_csv_file:
        csv_writer = csv.writer(new_csv_file, delimiter='\t')
        csv_writer.writerows(pair_ls)
    end = time.perf_counter()
    print(f'from {len(sorted_ls)} reads generate {len(pair_ls)} pairs, process takes {round(end - start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="Generate simulated pairs file from mapped barcode file from sprite_mapping.py")
    parser.add_argument("-in", help="mapped_barcode file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="simulated pairs file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()