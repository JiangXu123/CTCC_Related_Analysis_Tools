import csv
import seaborn as sns
import matplotlib.pyplot as plt


# to prepare the data for violin plot of the difference of the same_compartment_contact between CTCC and HiC
# plot was done in jupyter notebook under 'Matplotlib and Violin plot using Seaborn'
HiC_file = "/Users/jiangxu/Documents/Seq_data_analysis/dilution_HiC/merge_pair_out/hic_results/data/GM12878_hg19_MQ30/GM12878_diHiC_HindIII_hg19_MQ30_same_comp_contact_ratio.bed"
CTCC_file = "/Users/jiangxu/Documents/Seq_data_analysis/cmTCC/CTCC_SDS_HindIII_201908_MQ30_pair_compartment_ratio.bed"
merged_file = "/Users/jiangxu/Documents/Seq_data_analysis/cmTCC/CTCC_SDS_HindIII_HiC_HindIII_cis_trans_compare_for_plotting.tsv"
# HiC_df = pd.read_csv(HiC_file, delimiter='\t', names=['chr_name', 'start_pos', 'end_pos', 'sub_compartment', 'cis_S/N', 'trans_S/N', '1st_chr', '1st_chr_pos', '2nd_chr', '2nd_chr_pos', '3rd_chr', '3rd_chr_pos'])
# CTCC_df = pd.read_csv(HiC_file, delimiter='\t', names=['chr_name', 'start_pos', 'end_pos', 'sub_compartment', 'cis_S/N', 'trans_S/N', '1st_chr', '1st_chr_pos', '2nd_chr', '2nd_chr_pos', '3rd_chr', '3rd_chr_pos'])

with open(HiC_file, 'r') as hic_file:
    with open(CTCC_file, 'r') as ctcc_file:
        with open(merged_file, 'w') as m_file:
            hic_reader = csv.reader(hic_file, delimiter='\t')
            ctcc_reader = csv.reader(ctcc_file, delimiter='\t')
            m_writer = csv.writer(m_file, delimiter='\t')
            m_writer.writerow(['chr_name', 'start_pos', 'end_pos', 'sub_compartment', 'cis_sig_total', 'trans_sig_total', 'experiment'])
            while True:
                hic_line = next(hic_reader)
                ctcc_line = next(ctcc_reader)
                m_writer.writerow([hic_line[0], hic_line[1], hic_line[2], hic_line[3], hic_line[4], hic_line[5], 'diHiC'])
                m_writer.writerow([ctcc_line[0], ctcc_line[1], ctcc_line[2], ctcc_line[3], ctcc_line[4], ctcc_line[5], 'CTCC'])
