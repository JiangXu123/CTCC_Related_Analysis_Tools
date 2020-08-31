#! /usr/bin/env python

import argparse
import csv
import time


# to generate bin from  peak coordinate file
def run(args):
    start = time.perf_counter()
    bin_radius = args.radius
    peak_pos_file = args.ppf  # breaking position file: chr_name start end
    output_bin_file = args.output
    i = 0
    pre_chr = 'chr1'

    with open(peak_pos_file, 'r') as peak_pos:
        with open(output_bin_file, 'w') as bin_file:
            csv_writer = csv.writer(bin_file, delimiter='\t')
            csv_reader = csv.reader(peak_pos, delimiter='\t')
            pre_peak_pos = int(next(csv_reader)[2])
            for line in csv_reader:
                if pre_chr == line[0]:
                    if int(line[2]) >= pre_peak_pos + 2*bin_radius:  # to guarantee there's no bin overlap
                        pre_peak_pos = int(line[2])
                        left_boarder = (int(line[2]) - bin_radius)
                        right_boarder = (int(line[2]) + bin_radius)
                        csv_writer.writerow([line[0], left_boarder, right_boarder, float(line[3]), int(line[2]), str(i)])  # chr left_boarder right_boarder
                        i += 1  # chr_name, left boarder, right_boarder, height, center_pos,bin_number
                elif pre_chr != line[0]:  # if a new chr is encountered
                    pre_chr = line[0]
                    pre_peak_pos = int(line[2])

    end = time.perf_counter()
    print(f'{i} peak processed in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="to generate bin from  peak coordinate file")
    parser.add_argument("-r", help="bin radius", dest="radius", type=int, required=True)
    parser.add_argument("-ppf", help="peak position file", dest="ppf", type=str, required=True)
    parser.add_argument("-out", help="output bin file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()