import pandas as pd
import csv
import seaborn as sns
import matplotlib.pyplot as plt

HiC_file = "/Users/jiangxu/Documents/Seq_data_analysis/dilution_HiC/merge_pair_out/hic_results/data/GM12878_hg19_MQ30/GM12878_diHiC_HindIII_hg19_MQ30_same_comp_contact_ratio.bed"
CTCC_file = "/Users/jiangxu/Documents/Seq_data_analysis/cmTCC/CTCC_SDS_HindIII_201908_MQ30_pair_compartment_ratio.bed"
merged_file = "/Users/jiangxu/Documents/Seq_data_analysis/cmTCC/CTCC_SDS_HindIII_HiC_HindIII_hg19_10k.bed"
# HiC_df = pd.read_csv(HiC_file, delimiter='\t', names=['chr_name', 'start_pos', 'end_pos', 'sub_compartment', 'cis_S/N', 'trans_S/N', '1st_chr', '1st_chr_pos', '2nd_chr', '2nd_chr_pos', '3rd_chr', '3rd_chr_pos'])
# CTCC_df = pd.read_csv(HiC_file, delimiter='\t', names=['chr_name', 'start_pos', 'end_pos', 'sub_compartment', 'cis_S/N', 'trans_S/N', '1st_chr', '1st_chr_pos', '2nd_chr', '2nd_chr_pos', '3rd_chr', '3rd_chr_pos'])

with open(HiC_file, 'r') as hic_file:
    with open(CTCC_file, 'r') as ctcc_file:
        with open(merged_file, 'w') as m_file:
            hic_reader = csv.reader(hic_file, delimiter='\t')
            ctcc_reader = csv.reader(ctcc_file, delimiter='\t')
            m_writer = csv.writer(m_file, delimiter='\t')
            m_writer.writerow(['chr_name', 'start_pos', 'end_pos', 'sub_compartment', 'ctcc_cis', 'hic_cis', 'ctcc_trans', 'hic_trans', 'ctcc/hic(cis)', 'ctcc/hic(trans)'])
            while True:
                hic_line = next(hic_reader)
                ctcc_line = next(ctcc_reader)
                if (hic_line[0] == ctcc_line[0]) & (hic_line[1] == ctcc_line[1]) & (hic_line[2] == ctcc_line[2]):
                    ctcc_hic_cis_ratio = 0
                    ctcc_hic_trans_ratio = 0
                    hic_cis_value = float(hic_line[4])
                    hic_trans_value = float(hic_line[5])
                    ctcc_cis_value = float(ctcc_line[4])
                    ctcc_trans_value = float(ctcc_line[5])
                    if (hic_cis_value != 0) & (ctcc_cis_value != 0):
                        ctcc_hic_cis_ratio = round(ctcc_cis_value/hic_cis_value, 3)
                    if (hic_trans_value != 0) & (ctcc_trans_value != 0):
                        ctcc_hic_trans_ratio = round(ctcc_trans_value/hic_trans_value, 3)
                    compartment = ctcc_line[3]
                    m_writer.writerow([hic_line[0], hic_line[1], hic_line[2], compartment, ctcc_cis_value, hic_cis_value, ctcc_trans_value, hic_trans_value, ctcc_hic_cis_ratio, ctcc_hic_trans_ratio])







