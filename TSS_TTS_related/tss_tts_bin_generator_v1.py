import csv
import tabix

# To generate the bin for aggragate the possible interaction between TSS and TTS
radius_before_small = 20000
radius_after_high = 20000
feature = 'transcript'
chrom_size_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/hg19.chrom_ed.size'
compartment_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/HiC_comparmtent_file/GM12878_subcomp_sorted.bed.gz"
human_gff_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/gencode.v19.annotation.gff3.tsv"
TSS_bin_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/tss_tts_20k/gencode_v19_anno_gff3_tss_outer20k_bin.tsv"
TTS_bin_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/tss_tts_20k/gencode_v19_anno_gff3_tts_outer20k_bin.tsv"
feature_tss_tts_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/protein_coding_positive_strand_2000more_transcript_tss_tts_hg19.tsv'

bin_num = 1
chrom_dic = {}
chrom_ls = []

tb = tabix.open(compartment_file)
with open(chrom_size_file) as file1:
    csv_reader= csv.reader(file1, delimiter='\t')
    for line in csv_reader:
        chrom_dic.update({line[0]: int(line[1])})
        chrom_ls.append(line[0])
print(chrom_dic)
print(chrom_ls)

with open(human_gff_file, 'r') as file1:
    with open(TSS_bin_file, 'w') as file2:
        with open(TTS_bin_file, 'w') as file3:
            with open(feature_tss_tts_file, 'w') as file4:
                csv_reader = csv.reader(file1, delimiter='\t')
                csv_writer2 = csv.writer(file2, delimiter='\t')
                csv_writer3 = csv.writer(file3, delimiter='\t')
                csv_writer4 = csv.writer(file4, delimiter='\t')
                csv_writer4.writerow(['chr_name', 'TSS', 'TTS', 'length','strand', 'gene_type', 'compartment'])
                for line in csv_reader: # the original annotation file
                    if '##sequence' in line[0]:  # these are the leading title for chromosomes
                        pass
                    elif line[2] == feature:  # look for feature = 'transcript'
                        try:
                            low = int(line[3])
                            high = int(line[4])
                            chr_name = line[0]
                            low_query_result = tb.query(chr_name, low-1, low)
                            high_query_result = tb.query(chr_name, high-1, high)
                            low_compartment = next(low_query_result)[3][:2]
                            high_compartment = next(high_query_result)[3][:2]  # get the compartment belong to that TSS_TTS
                            strand = line[6]
                            gene_type = line[8].split(';')[4].split('=')[1]
                            transcript_length = high - low
                            if (gene_type == 'protein_coding') & (high_compartment == 'A1') & (transcript_length > 2000) & (strand == '+'):
                                csv_writer4.writerow([line[0], int(line[3]), int(line[4]), transcript_length, strand, gene_type, low_compartment, high_compartment])

                            inner_diameter = int((high - low)/2)
                            high_upper_boarder = high + radius_after_high
                            high_inner_boarder = high - inner_diameter
                            low_lower_boarder = low - radius_before_small
                            low_inner_boarder = low + inner_diameter
                            if (high_upper_boarder < chrom_dic[line[0]]) & (low_lower_boarder > 0):
                                csv_writer2.writerow([chr_name, low_lower_boarder, low_inner_boarder, low, strand, gene_type, bin_num, low_compartment])
                                csv_writer3.writerow([chr_name, high_inner_boarder, high_upper_boarder, high, strand, gene_type, bin_num, high_compartment])
                                bin_num += 1
                            else:
                                pass
                        except:
                            pass
