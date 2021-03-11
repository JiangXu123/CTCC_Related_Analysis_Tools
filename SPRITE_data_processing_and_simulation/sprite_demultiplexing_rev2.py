#! /usr/bin/env python

import pandas as pd
import argparse
import time
import re
import gzip
# This code is used to extract the barcode information from the pairwise Sequencing file R1_fastq.gz and R2_fastq.gz
# The code has been optimized for running time, i.e. I didn't use the time-consuming pandas df.concat for data accumulation
# Instead, I used the ls.append() method. which is much faster.


def run(args):
    start = time.perf_counter()
    seq_R1_file = args.R1  # seq_file_R1 contains the DPM barcode and the genomic DNA sequencing information
    seq_R2_file = args.R2  # seq_file_R2 contains the odd, even and terminal barcode information
    output_file = args.output  # the barcode-tagged reads number tsv file
    DPM_bc_file = args.DPM  # the DPM barcode file
    Ter_bc_file = args.ter  # the Terminal barcode file
    Odd_bc_file = args.odd  # the Odd barcode file
    Even_bc_file = args.even  # the Even barcode file

    dpm_ls = []
    ter_ls = []
    odd_ls = []
    even_ls = []
    reads_info_ls = []
    DPM_bc_count = 0
    ter_bc_count = 0
    odd_bc_2_count = 0
    even_bc_count = 0
    odd_bc_1_count = 0
    all_bc_read_count = 0

    with open(DPM_bc_file, 'r') as dpm:
        i = 1
        for line in dpm.read().splitlines():
            dpm_ls.append([line, f'{i:02}'])
            i += 1
    # print(dpm_ls)

    with open(Ter_bc_file, 'r') as ter:
        i = 1
        for line in ter.read().splitlines():
            ter_ls.append([line, f'{i:02}'])
            i += 1
            # print(ter_ls)

    with open(Odd_bc_file, 'r') as odd:
        i = 1
        for line in odd.read().splitlines():
            odd_ls.append([line, f'{i:02}'])
            i += 1
            # print(odd_ls)

    with open(Even_bc_file, 'r') as even:
        i = 1
        for line in even.read().splitlines():
            even_ls.append([line, f'{i:02}'])
            i += 1
            # print(even_ls)

    with gzip.open(seq_R1_file, 'rt') as f1:  # open a gzipped file in text mode
        with gzip.open(seq_R2_file, 'rt') as f2:
            lines_1 = f1.read().splitlines()  # generate a list from the sequencing file R1
            lines_2 = f2.read().splitlines()  # generate a list from the sequencing file R2
            total_count = len(lines_1)/4
            print(total_count)
            for i in range(0, int(len(lines_1) / 4)):
                DPM_bc = f'{0:02}'
                odd_bc_1 = f'{0:02}'
                even_bc = f'{0:02}'
                odd_bc_2 = f'{0:02}'
                ter_bc = f'{0:02}'
                read_ls = []
                read_ls.append(lines_1[4 * i])  # get the reads number
                for barcode in dpm_ls:
                    pattern_DPM = re.compile(barcode[0])
                    if pattern_DPM.match(lines_1[4 * i + 1]):
                        DPM_bc_count += 1
                        DPM_bc = barcode[1]
                read_ls.append(DPM_bc)

                for barcode in ter_ls:
                    pattern_terminal = re.compile(barcode[0] + 'GACAACT')
                    if pattern_terminal.match(lines_2[4 * i + 1]):
                        ter_bc_count += 1
                        ter_bc = barcode[1]
                read_ls.append(ter_bc)

                for barcode in odd_ls:
                    pattern_sec_odd = re.compile(r'\w{16,19}' + barcode[0])
                    pattern_first_odd = re.compile(r'\w{65,68}' + barcode[0])
                    if pattern_sec_odd.match(lines_2[4 * i + 1]):
                        odd_bc_2_count += 1
                        odd_bc_2 = barcode[1]
                    if pattern_first_odd.match(lines_2[4 * i + 1]):
                        odd_bc_1_count += 1
                        odd_bc_1 = barcode[1]
                read_ls.append(odd_bc_2)
                read_ls.append(odd_bc_1)

                for barcode in even_ls:
                    pattern_even = re.compile(r'\w{41,44}' + barcode[0])
                    if pattern_even.match(lines_2[4 * i + 1]):
                        even_bc_count += 1
                        even_bc = barcode[1]
                read_ls.append(even_bc)

                if (read_ls[1] != '00') & (read_ls[2] != '00') & (read_ls[3] != '00') & (read_ls[4] != '00') & (
                        read_ls[5] != '00'):
                    all_bc_read_count += 1

                read_ls.append(''.join(read_ls[1:6]))
                reads_info_ls.append(read_ls)
    df = pd.DataFrame(reads_info_ls,
                      columns=['read_name', 'DPM_barcode', 'terminal_barcode', 'sec_odd_barcode', 'first_odd_barcode',
                               'even_barcode', 'complex_barcode'])
    print(f'total count is {total_count}')
    print(f'{100 * DPM_bc_count / total_count} percent reads have DPM barcode')
    print(f'{100 * ter_bc_count / total_count} percent reads have Ter barcode')
    print(f'{100 * odd_bc_2_count / total_count} percent reads have second odd barcode')
    print(f'{100 * even_bc_count / total_count} percent reads have even barcode')
    print(f'{100 * odd_bc_1_count / total_count} percent reads have first odd barcode')
    print(f'{100 * all_bc_read_count / total_count} percent reads have all 5 barcodes detected')
    end = time.perf_counter()
    print(f'{total_count} reads takes {round(end - start, 2)} seconds')
    df.to_csv(output_file, header=False, index=False, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="Demultiplexing SPRITE barcode")
    parser.add_argument("-R1", help="input R1 sequencing file", dest="R1", type=str, required=True)
    parser.add_argument("-R2", help="input R2 sequencing file", dest="R2", type=str, required=True)
    parser.add_argument("-DPM", help="DPM barcode file", dest="DPM", type=str, required=True)
    parser.add_argument("-ter", help="Terminal barcode file", dest="ter", type=str, required=True)
    parser.add_argument("-odd", help="Odd barcode file", dest="odd", type=str, required=True)
    parser.add_argument("-even", help="Even barcode file", dest="even", type=str, required=True)
    parser.add_argument("-out", help="the barcodes file generated", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()