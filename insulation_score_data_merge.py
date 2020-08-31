import csv

xlink_first_CTCC_SDS = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/insulation_score_correlation_analyis/xlink_first_CTCC_SDS_hg19_MQ0_10k_bl_100pixel_insulation.tsv"
xlink_second_CTCC_SDS = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/insulation_score_correlation_analyis/xlink_second_CTCC_SDS_27547k_10k_bl_100pixel_insulation.tsv"
output_file = '/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/insulation_score_correlation_analyis/xlink_first_second_insulation_10k_bl.tsv'

with open(xlink_first_CTCC_SDS, 'r') as file1:
    with open(xlink_second_CTCC_SDS, 'r') as file2:
        csv_reader1 = csv.reader(file1, delimiter='\t')
        csv_reader2 = csv.reader(file2, delimiter='\t')
        next(csv_reader1)
        next(csv_reader2)
        with open(output_file, 'w') as file3:
            csv_writer = csv.writer(file3, delimiter='\t')
            csv_writer.writerow(['chrom', 'xlink_first_insulation_score', 'xlink_second_insulation_score'])
            for line1, line2 in zip(csv_reader1, csv_reader2):
                try:
                    if (line1[0] == line2[0]) & (line1[1] == line2[1]):
                        csv_writer.writerow([line1[0], line1[4], line2[4]])
                        print(['chrom',line1[4], line2[4]])
                except IndexError:
                    pass

                # print(line2[0])

                # if (line1[0] == line2[0]) & (line1[1] == line2[1]): # if the bins from both file are the same
                #     csv_writer.writerow(['chrom', bin_number, line1[5], line2[5]])
                #     bin_number += 1





