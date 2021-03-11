import csv
import tabix

# To generate the bin for aggregating the possible interaction between TSS and TTS
bin_radius = 300
feature = 'transcript'
chrom_size_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/hg19.chrom_ed.size'
compartment_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/HiC_comparmtent_file/GM12878_subcomp_sorted.bed.gz"
human_gff_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/gencode.v19.annotation.gff3.tsv"
TSS_bin_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/gencode_v19_anno_gff3_tss_r300_bin.tsv"
TTS_bin_file = "/Users/jiangxu/Documents/Seq_data_analysis/genome_file/gencode_v19_anno_gff3_tts_r300_bin.tsv"
feature_tss_tts_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/transcript_tss_tts_hg19.tsv'

bin_num = 1
chrom_dic = {}
chrom_ls = []

tb = tabix.open(compartment_file)
with open(chrom_size_file) as file1:
    csv_reader = csv.reader(file1, delimiter='\t')
    for line in csv_reader:
        chrom_dic.update({line[0]:int(line[1])})
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
                csv_writer4.writerow(['chr_name', 'TSS', 'TTS', 'length', 'strand', 'gene_type', 'compartment'])
                for line in csv_reader: # the original annotation file
                    if '##sequence' in line[0]:  # these are the leading title for chromosomes
                        pass
                    elif line[2] == feature:  # look for feature = 'transcript'
                        try:
                            query_result = tb.query(line[0], int(line[3]), int(line[4]))
                            compartment = next(query_result)[3]  # get the compartment belong to that TSS_TTS
                            strand = line[6]
                            gene_type = line[8].split(';')[4].split('=')[1]
                            transcript_length = int(line[4]) - int(line[3])
                            csv_writer4.writerow([line[0], line[3], line[4], transcript_length, strand, gene_type, compartment])
                            if strand == '+':  # if the transcript is from the forward strand
                                tss = int(line[3])
                                tts = int(line[4])
                                if tts + bin_radius < chrom_dic[line[0]]:
                                    tss_bin_start = tss - bin_radius
                                    tss_bin_end = tss + bin_radius
                                    tts_bin_start = tts - bin_radius
                                    tts_bin_end = tts + bin_radius
                                elif tts + bin_radius >= chrom_dic[line[0]]:
                                    if tss + bin_radius < chrom_dic[line[0]]:
                                        tss_bin_start = tss - bin_radius
                                        tss_bin_end = tss + bin_radius
                                        tts_bin_start = tts - bin_radius
                                        tts_bin_end = chrom_dic[line[0]]
                                    elif tss + bin_radius >= chrom_dic[line[0]]:
                                        tss_bin_start = tss - bin_radius
                                        tss_bin_end = chrom_dic[line[0]]
                                        tts_bin_start = tts - bin_radius
                                        tts_bin_end = chrom_dic[line[0]]
                                csv_writer2.writerow([line[0], tss_bin_start, tss_bin_end, tss, strand, gene_type, bin_num, compartment])
                                csv_writer3.writerow([line[0], tts_bin_start, tts_bin_end, tts, strand, gene_type, bin_num, compartment])
                                bin_num += 1
                            elif strand == '-':  # if the transcript is from the reverse strand
                                tss = int(line[4])
                                tts = int(line[3])
                                if tss + bin_radius < chrom_dic[line[0]]:
                                    tss_bin_start = tss + bin_radius
                                    tss_bin_end = tss - bin_radius
                                    tts_bin_start = tts + bin_radius
                                    tts_bin_end = tts - bin_radius
                                elif tss + bin_radius >= chrom_dic[line[0]]:
                                    if tts + bin_radius < chrom_dic[line[0]]:
                                        tss_bin_start = chrom_dic[line[0]]
                                        tss_bin_end = tss - bin_radius
                                        tts_bin_start = tts + bin_radius
                                        tts_bin_end = tts - bin_radius
                                    elif tts + bin_radius >= chrom_dic[line[0]]:
                                        tss_bin_start = chrom_dic[line[0]]
                                        tss_bin_end = tss - bin_radius
                                        tts_bin_start = chrom_dic[line[0]]
                                        tts_bin_end = tts - bin_radius
                                csv_writer2.writerow([line[0], tts_bin_end, tts_bin_start, tts, strand, gene_type, bin_num, compartment])
                                csv_writer3.writerow([line[0], tss_bin_end, tss_bin_start, tss, strand, gene_type, bin_num, compartment])
                                bin_num += 1
                        except:
                            pass
