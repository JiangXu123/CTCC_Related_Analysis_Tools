#! /usr/bin/env python

import argparse
import csv
import pandas as pd
import time


def run(args):
    pair_dist_freq_file = args.input
    output_pair_compartment_file = args.output
    exp_1_name = args.name_1
    exp_2_name = args.name_2
    df = pd.read_csv(pair_dist_freq_file, delimiter='\t')
    df1 = df.loc[df['Experiment'] == exp_1_name]
    df2 = df.loc[df['Experiment'] == exp_2_name]
    df1.to_csv(f"{exp_1_name}.csv", sep='\t')
    df2.to_csv(f"{exp_2_name}.csv", sep='\t')
    accumulated_freq_dif = 0
    print(f"DataFrame files {exp_1_name}.csv and {exp_1_name}.csv generated")
    with open(f"{exp_1_name}.csv", "r") as file1:
        with open(f"{exp_2_name}.csv", "r") as file2:
            with open(output_pair_compartment_file, "w") as file3:
                csv_reader1 = csv.reader(file1, delimiter='\t')
                csv_reader2 = csv.reader(file2, delimiter='\t')
                csv_writer = csv.writer(file3, delimiter='\t')
                next(csv_reader1)
                next(csv_reader2)
                line1 = next(csv_reader1)
                line2 = next(csv_reader2)
                csv_writer.writerow(["Distance", "Difference_in_Counts_Per_Million", "Accumulated_Difference_in_Counts_Per_Million"])
                while True:
                    try:
                        while int(line1[1]) > int(line2[1]):
                            line2 = next(csv_reader2)
                        while int(line1[1]) < int(line2[1]):
                            line1 = next(csv_reader1)
                        if line1[1] == line2[1]:  # if the distance is the same
                            difference = int(line1[3]) - int(line2[3])  # the difference in counts per million
                            accumulated_freq_dif += difference
                            csv_writer.writerow([line1[1], difference, accumulated_freq_dif])
                            line1 = next(csv_reader1)
                            line2 = next(csv_reader2)
                    except:
                        pass
                    # print(f"line1's distance is {line1[1]}and line2's distance is {line2[1]}")
    end = time.perf_counter()
    print(f'process finished in {end-start} sec')


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-i", help="the input cis_dist_file", dest="input", type=str, required=True)
    parser.add_argument("-e1", help="the name of first experiment", dest="name_1", type=str, required=True)
    parser.add_argument("-e2", help="the name of second experiment", dest="name_2", type=str, required=True)
    parser.add_argument("-o", help="output files name", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
