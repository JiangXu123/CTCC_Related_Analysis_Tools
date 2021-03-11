#! /usr/bin/env python

import argparse
import csv
import time


# calculate the center location of peaks from a bedgraph file
# generate the bedgraph file from bigwigfile using bigWigToBedGraph
def run(args):
    start = time.perf_counter()
    input_file = args.input
    output_file = args.output
    high_threshold = args.th
    low_threshold = args.tl
    feature = args.feature
    gp = args.gap  # the default value, individual peak who are close(but not continuous) will not be merged
    temp_ls = []  # used for store only peak related lines
    pre_slop = 0.5  # pre_slop is arbitrarily set to 0.5
    with open(input_file, 'r') as file1:
        with open(output_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            temp_ls.append(next(csv_reader))  # put the the first line into temp list
            csv_writer = csv.writer(file2, delimiter='\t')
            for line in csv_reader:  # parsing start from the second line
                previous_line_chr = temp_ls[-1][0]
                previous_line_pos = int(temp_ls[-1][2])
                previous_line_value = float(temp_ls[-1][3])
                new_line_chr = line[0]
                new_line_pos = int(line[1])
                new_line_value = float(line[3])

                if new_line_pos - previous_line_pos > gp:  # if new_line belong to another separate peak
                    pre_slop = 0.5
                    temp_ls = [line]
                elif previous_line_chr != new_line_chr:  # if new_line belong to another chromosome
                    pre_slop = 0.5
                    print(f'{previous_line_chr} finished')
                    temp_ls = [line]

                elif ((new_line_pos - previous_line_pos) == gp) & (previous_line_chr == new_line_chr):
                    # if it is within the same chromosome and is continuous with previous position
                    if (new_line_value > previous_line_value) & (pre_slop > 0):
                        temp_ls = [line]  # this is the feature of an up slope. temp_ls is emptied and present line is appended.
                        pre_slop = new_line_value - previous_line_value

                    elif (new_line_value < previous_line_value) & (pre_slop > 0) & (len(temp_ls) == 1):  # if it is a peak top
                        peak_width = (int(temp_ls[0][2]) - int(temp_ls[0][1]))
                        center_pos = (int(temp_ls[0][1]) + round(peak_width / 2))  # if peak_width is 1, then round(1/2)=0
                        # so both flat peak and sharp peak are considered
                        pre_slop = new_line_value - previous_line_value
                        if (feature == 'peak') & (previous_line_value >= low_threshold) & (previous_line_value <= high_threshold):
                            csv_writer.writerow([temp_ls[0][0], center_pos - 1, center_pos, previous_line_value])
                        temp_ls = [line]

                    elif (new_line_value < previous_line_value) & (pre_slop < 0) & (len(temp_ls) == 1):
                        pre_slop = new_line_value - previous_line_value  # now the curve continue to drop
                        temp_ls = [line]

                    elif (new_line_value > previous_line_value) & (pre_slop < 0) & (len(temp_ls) == 1):
                        peak_width = (int(temp_ls[0][2]) - int(temp_ls[0][1]))
                        center_pos = (int(temp_ls[0][1]) + round(peak_width / 2))  # this is a concave flat valley
                        pre_slop = new_line_value - previous_line_value  # now the curve begin to rise
                        if (feature == 'valley') & (previous_line_value >= low_threshold) & (previous_line_value <= high_threshold):  # if you want to write the valley position to a file
                            csv_writer.writerow([temp_ls[0][0], center_pos - 1, center_pos, previous_line_value])
                        temp_ls = [line]
    end = time.perf_counter()
    print(f'process takes {round(end-start, 2)} seconds')


def main():
    parser = argparse.ArgumentParser(description="calculate the center location of peaks or valleys in bedgraph file")
    parser.add_argument("-i", help="input bedgraph file)", dest="input", type=str, required=True)
    parser.add_argument("-ft", help="feature to calculate, peak or valley?", dest="feature", type=str, required=True)
    parser.add_argument("-th", help="high value threshold", dest="th", type=float, required=True)
    parser.add_argument("-tl", help="low value threshold", dest="tl", type=float, required=True)
    parser.add_argument("-g", help="gap value, separate peaks within the gap are regarded as a single peak", dest="gap", type=int, required=True)
    parser.add_argument("-o", help="peak position file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()