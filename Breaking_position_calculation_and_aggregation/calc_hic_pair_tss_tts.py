#! /usr/bin/env python

import argparse
import tabix
import csv
import time


# aggregate breaking position relative to certain tabix indexed peak position file,
def run(args):
    pair_file = args.pf
    tss_bin_file = args.tssb
    tts_bin_file = args.ttsb
    tss_tts_aggr_file = args.output

    tb_tss = tabix.open(tss_bin_file)  # tss_bin_file and tts_bin_file are both sorted bgzip compressed
    tb_tts = tabix.open(tts_bin_file)  # and tabix indexed file
    line_number = 0

    start = time.perf_counter()
    with open(pair_file, 'r') as file1:
        with open(tss_tts_aggr_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            for line in csv_reader:
                r1_chr = line[1]
                r1_pos = int(line[2])
                r2_chr = line[4]
                r2_pos = int(line[5])
                try:
                    tss_results = tb_tss.query(r1_chr, r1_pos - 1, r1_pos)  # in .ValidPair file,
                    tts_results = tb_tts.query(r2_chr, r2_pos - 1, r2_pos)
                    for tss_result in tss_results:
                        for tts_result in tts_results:
                            if tss_result[6] == tts_result[6]:  # if a reads pair fall into the same TSS_TTS pair, and in the same chromosome.
                                compartment_tss = tss_result[7][:2]
                                bin_num = int(tss_result[6])
                                bin_sign = tss_result[4]
                                gene_type = tss_result[5]
                                compartment_tts = tts_result[7][:2]
                                tss_cen_pos = int(tss_result[3])
                                tts_cen_pos = int(tts_result[3])
                                read1_to_tss_cen_distance = r1_pos - tss_cen_pos
                                read2_to_tts_cen_distance = r2_pos - tts_cen_pos
                                feature_length = tts_cen_pos - tss_cen_pos
                                if bin_sign == '+':  # if the reference bin is on the forward strand
                                    csv_writer.writerow([line[1], read1_to_tss_cen_distance, read2_to_tts_cen_distance, bin_sign, bin_num, gene_type, feature_length, compartment_tss, compartment_tts])
                                if bin_sign == '-':  # if the reference bin is on the reverse strand
                                    csv_writer.writerow([line[1], -read2_to_tts_cen_distance, -read1_to_tss_cen_distance, bin_sign, bin_num, gene_type, feature_length, compartment_tts, compartment_tss])
                        # if reference bin is on the. reverse strand, because read2 is always mapped to higher position in .validPairs from HiC-Pro
                        # than read1, read2_to_tts_cen_distance become the distance_to_tss on reverse strand, read1_to_tss_cen_distance become the distance_to_tts on reverse strand
                        # also, a positive distance to the center(either TSS or TTS) always means a negative distance to the center on the reverse strand
                        # The columns are: chr_name, dist_to_tss, dist_to_tts, bin_sign, bin_num, gene_type, compartment_tss, compartment_tts
                except:
                    pass
                line_number += 1
                if line_number % 10000 == 0:
                    print(f'{line_number} of pairs processed')

    end = time.perf_counter()
    print(f'{line_number} breaking position processed in {round(end - start, 2)} sec')


def main():
    parser = argparse.ArgumentParser(description="to calculate HiC pairs distribution relative to a pannel of TSS and TTS")
    parser.add_argument("-tssb", help="bgzipped tabix indexed TSS bin file", dest="tssb", type=str, required=True)
    parser.add_argument("-ttsb", help="bgzipped tabix indexed TTS bin file", dest="ttsb", type=str, required=True)
    parser.add_argument("-pf", help="HiC-Pro .allValidPair file", dest="pf", type=str, required=True)
    parser.add_argument("-out", help="calculated list file", dest="output", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()