import numpy as np
import csv
import matplotlib.pyplot as plt

input_matrix1 = '/Users/jiangxu/Documents/Seq_data_analysis/tutorials/coolpuppy_tutorial/Ctrl_loop.tsv'
input_matrix2 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/with_with_out_RNaseA/CTCC_SDS_noRNaseA_size_matched_hg19_200b_CTCF_related_header_removed.tsv'
input_matrix3 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/with_with_out_RNaseA/CTCC_SDS_noRNaseA_size_matched_hg19_200b_50peak_CTCF_related_header_removed.txt'
input_matrix4 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/with_with_out_RNaseA/CTCC_SDS_noRNaseA_size_matched_hg19_200b_50peak_1000_pad_CTCF_related_header_removed.txt'
input_matrix5 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/with_with_out_RNaseA/CTCC_SDS_noRNaseA_size_matched_hg19_200b_50peak_100_pad_CTCF_related_header_removed.txt'
input_matrix6 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/with_with_out_RNaseA/test1_header_removed.txt'
input_matrix7 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/GM12878_SDS_noRNaseA_CTCF_pileup_pad_250_header_removed.tsv'
input_matrix8 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/GM12878_SDS_withRNaseA_CTCF_pileup_pad_250_header_removed.tsv'
input_matrix9 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/CTCC_all_hg19_MQ0_100bp_TSS_TTS_plustrand_apileup_header_romved.tsv'
input_matrix10 = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/CTCC_all_hg19_MQ0_100bp_TSS_TTS_plustrand_apileup_header_removed.tsv'
CTCC1_matrix = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/with_without_cryomilling/CTCC1_SDS_RaoPrimary_loop_pileup_100bp.matrix'
H = []
line_ls = []
with open(CTCC1_matrix, 'r') as file1:
    csv_reader = csv.reader(file1, delimiter=' ')
    for line in csv_reader:
        print(line)
        for item in line:
            line_ls.append(float(item))
        H.append(line_ls)
        line_ls = []
print(H)
# H = np.array([[1, 2, 3, 4],
#               [5, 6, 7, 8],
#               [9, 10, 11, 12],
#               [13, 14, 15, 16]])  # added some commas and array creation code

fig = plt.figure(figsize=(6, 3.2))

ax = fig.add_subplot(111)
ax.set_title('CTCF peak centered pile up')
plt.imshow(H)
ax.set_aspect('equal')

cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()