#! /usr/bin/env python

import argparse
import tabix
import csv


input_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/mouse/gencode.vM25.annotation.gff3.tsv'
# compartment_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/HiC_comparmtent_file/GM12878_subcomp_sorted.bed.gz'
tss_tts_bed_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/mouse/mm10_tss_tts_transcript_2klength_plusstrand.tsv'
length_lower_th = 2000
length_higher_th = 20000
# compartment = 'A1'   # possible choice: A1, B1, B2, B3, B4, NA
feature_strand = '+'
feature = 'transcript'
# tb = tabix.open(compartment_file)
previous_high = 0

with open(input_file, 'r') as file1:
    with open(tss_tts_bed_file, 'w') as file2:
        csv_reader = csv.reader(file1, delimiter='\t')
        csv_writer2 = csv.writer(file2, delimiter='\t')
        for line in csv_reader: # the transcript annotation file
            if '##sequence' not in line[0]:
                if line[2] == feature:
                    feature_length = int(line[4]) - int(line[3])
                    strand = line[6]
                    low = int(line[3])
                    high = int(line[4])
                    chr_name = line[0]
                    # try:
                    #     low_query_result = tb.query(chr_name, low - 1, low)
                    #     low_compartment = next(low_query_result)[3][:2]  # get the compartment belong to that TSS_TTS
                    if (feature_length >= length_lower_th) & (feature_length <= length_higher_th) & (feature_strand == strand) & (low > previous_high):
                        print(line)
                        csv_writer2.writerow([chr_name, low, high])
                    previous_high = high

