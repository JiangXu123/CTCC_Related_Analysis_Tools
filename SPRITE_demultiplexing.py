#! /usr/bin/env python

import pandas as pd
import argparse
import time
import re


def run(args):
    start = time.perf_counter()
    seq_R1_file = args.R1  # seq_file_R1 contains the DPM barcode and the genomic DNA sequencing information
    seq_R2_file = args.R2  # seq_file_R2 contains the odd, even and terminal barcode information
    output_file = args.output  # the barcode-tagged reads number tsv file
    DPM_bc_file = args.DPM  # the DPM barcode file
    Ter_bc_file = args.ter  # the Terminal barcode file
    Odd_bc_file = args.odd  # the Odd barcode file
    Even_bc_file = args.even  # the Even barcode file

    seq_R1_df = pd.read_csv(seq_R1_file, header=None)  # seq_R1_df contains the odd, even, odd, terminal barcode information
    seq_R2_df = pd.read_csv(seq_R2_file, header=None)  # seq_R2_df contains the DPM barcoding information and

    DPM_bc_df = pd.read_csv(DPM_bc_file, header=None)
    DPM_bc_ls = []
    i = 1
    for index, row in DPM_bc_df.iterrows():
        DPM_bc_ls.append([f'{i:02}', row[0]])
        i += 1
    print(DPM_bc_ls)

    ter_bc_df = pd.read_csv(Ter_bc_file, names=['barcode_name', 'barcode'], delimiter='\t')
    ter_bc_ls = []
    i = 1
    for index, row in ter_bc_df.iterrows():
        ter_bc_ls.append([f'{i:02}', row[1]])
        i += 1
    print(ter_bc_ls)

    odd_bc_df = pd.read_csv(Odd_bc_file, names=['barcode'], delimiter='\t')
    odd_bc_ls = []
    i = 1
    for index, row in odd_bc_df.iterrows():
        odd_bc_ls.append([f'{i:02}', row[0]])
        i += 1
    print(odd_bc_ls)

    even_bc_df = pd.read_csv(Even_bc_file, names=['barcode'], delimiter='\t')
    even_bc_ls = []
    i = 1
    for index, row in even_bc_df.iterrows():
        even_bc_ls.append([f'{i:02}', row[0]])
        i += 1
    print(even_bc_ls)

    DPM_bc_count = 0
    ter_bc_count = 0
    odd_bc_2_count = 0
    even_bc_count = 0
    odd_bc_1_count = 0
    non_DPM_bc_count = 0
    non_ter_bc_count = 0
    non_odd_bc_2_count = 0
    non_even_bc_count = 0
    non_odd_bc_1_count = 0
    all_bc_read_count = 0
    reads_info_ls = []
    print(len(seq_R2_df) / 4)
    for i in range(0, int(len(seq_R2_df) / 4)):
        DPMbc_unmatch_i = 0
        terbc_unmatch_i = 0
        second_oddbc_unmatch_i = 0
        evenbc_unmatch_i = 0
        first_oddbc_unmatch_i = 0
        DPM_bc = f'{0:02}'
        odd_bc_1 = f'{0:02}'
        even_bc = f'{0:02}'
        odd_bc_2 = f'{0:02}'
        ter_bc = f'{0:02}'
        read_ls = []
        read_ls.append(seq_R1_df.iloc[4 * i, 0])  # get the reads number
        for barcode in DPM_bc_ls:
            pattern_DPM = re.compile(barcode[1])
            if pattern_DPM.match(seq_R1_df.iloc[4 * i + 1, 0]):
                DPM_bc_count += 1
                DPM_bc = barcode[0]
            else:
                DPMbc_unmatch_i += 1
                if DPMbc_unmatch_i == 96:
                    non_DPM_bc_count += 1
        read_ls.append(DPM_bc)

        for barcode in ter_bc_ls:
            pattern_terminal = re.compile(barcode[1] + 'GACAACT')
            if pattern_terminal.match(seq_R2_df.iloc[4 * i + 1, 0]):
                ter_bc_count += 1
                ter_bc = barcode[0]
            else:
                terbc_unmatch_i += 1
                if terbc_unmatch_i == 96:
                    non_ter_bc_count += 1
        read_ls.append(ter_bc)

        for barcode in odd_bc_ls:
            pattern_sec_odd = re.compile(r'\w{16,19}' + barcode[1])
            if pattern_sec_odd.match(seq_R2_df.iloc[4 * i + 1, 0]):
                odd_bc_2_count += 1
                odd_bc_2 = barcode[0]
            else:
                second_oddbc_unmatch_i += 1
                if second_oddbc_unmatch_i == 96:
                    non_odd_bc_2_count += 1
        read_ls.append(odd_bc_2)

        for barcode in even_bc_ls:
            pattern_even = re.compile(r'\w{41,44}' + barcode[1])
            if pattern_even.match(seq_R2_df.iloc[4 * i + 1, 0]):
                even_bc_count += 1
                even_bc = barcode[0]
            else:
                evenbc_unmatch_i += 1
                if evenbc_unmatch_i == 96:
                    non_even_bc_count += 1
        read_ls.append(even_bc)

        for barcode in odd_bc_ls:
            pattern_first_odd = re.compile(r'\w{65,68}' + barcode[1])
            if pattern_first_odd.match(seq_R2_df.iloc[4 * i + 1, 0]):
                odd_bc_1_count += 1
                odd_bc_1 = barcode[0]
            else:
                first_oddbc_unmatch_i += 1
                if first_oddbc_unmatch_i == 96:
                    non_odd_bc_1_count += 1
        read_ls.append(odd_bc_1)
        if (read_ls[1] != '00') & (read_ls[2] != '00') & (read_ls[3] != '00') & (read_ls[4] != '00') & (read_ls[5] != '00'):
            all_bc_read_count += 1
            new_ls = [read_ls[0], ''.join(read_ls[1:6])]
            reads_info_ls.append(new_ls)
    reads_info_ls.append(new_ls)
    output_df = pd.DataFrame(reads_info_ls, columns=['read_name', 'complex_barcode'])
    total_count = DPM_bc_count + non_DPM_bc_count
    print(f'{100 * DPM_bc_count / total_count} percent reads have DPM barcode')
    print(f'{100 * ter_bc_count / total_count} percent reads have Ter barcode')
    print(f'{100 * odd_bc_2_count / total_count} percent reads have second odd barcode')
    print(f'{100 * even_bc_count / total_count} percent reads have even barcode')
    print(f'{100 * odd_bc_1_count / total_count} percent reads have first odd barcode')
    print(f'{100 * all_bc_read_count / total_count} percent reads have all 5 barcodes detected')
    end = time.perf_counter()
    print(f'{total_count} reads takes {round(end - start, 2)} seconds')

    output_df.to_csv(output_file, sep='\t', index=False, header=None)


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