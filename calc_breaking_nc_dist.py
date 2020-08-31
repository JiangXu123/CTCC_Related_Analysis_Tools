#! /usr/bin/env python

import argparse
import csv
import time
import tabix


# extract breaking information from the breaking position with compartment file
# and calculate the distance to the center of the nucleosome
def run(args):
    start = time.perf_counter()
    peak_sig_indexed_bin_file = args.pscib
    bin_radius = args.radius
    break_pos_file = args.bpf  # breaking position file: chr_name start end
    output_aggr_ls = args.output
    exp = args.exp
    tb = tabix.open(peak_sig_indexed_bin_file)
    line_count = 0
    with open(break_pos_file, 'r') as break_pos:
        with open(output_aggr_ls, 'w') as aggregated_data:
            csv_writer = csv.writer(aggregated_data, delimiter='\t')
            csv_reader = csv.reader(break_pos, delimiter='\t')
            for line in csv_reader:
                try:
                    result = tb.query(line[0], int(line[1]), int(line[2])) # result1 contains: chr_name, left boarder, right_boarder, height, center_pos,bin_number
                    item = next(result)
                    distance = int(line[2]) - int(item[4])  # item[4] contains the bin center position
                    if abs(distance) <= bin_radius:
                        csv_writer.writerow([exp, int(item[5]), line[3], line[0], distance])
                        # columns name are: experiment name, nucleosome number(bin number), sub_nucleus compartment, chr_name,
                        # distance to nucleosome center
                except:
                    pass
                line_count += 1
                if line_count % 100000 == 0:
                    print(f'{line_count} of lines processed')

    end = time.perf_counter()
    print(f'{line_count} breaking position processed in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="to calculate breaking frequecy relative to a buch of peaks center position")
    parser.add_argument("-pscib", help="peak signal center indexed bin file", dest="pscib", type=str, required=True)
    parser.add_argument("-r", help="bin radius", dest="radius", type=int, required=True)
    parser.add_argument("-bpf", help="breaking position file", dest="bpf", type=str, required=True)
    parser.add_argument("-expt", help="experiment name", dest="exp", type=str, required=True)
    parser.add_argument("-out", help="aggregated list file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()