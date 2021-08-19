import csv

bin_size = 10000
chromosome_ls = []
chromosome_size_file = '/Users/jiangxu/Documents/Seq_data_analysis/genome_file/hg19.chrom_ed.size'
comp_tagged_pair_file = '/Users/jiangxu/Documents/Seq_data_analysis/dilution_HiC/merge_pair_out/hic_results/data/GM12878_hg19_MQ30/GM12878_diHiC_HindIII_hg19_MQ30_comp_pair_tagged.tsv'
chromosome_size_dic = {}
chr_info_dic = {}
chr_info_sorted_ls = []
with open(chromosome_size_file, 'r') as chr_size_file:
    csv_reader = csv.reader(chr_size_file, delimiter='\t')
    for line in csv_reader:
        chromosome_ls.append(line[0])
        chromosome_size_dic[line[0]] = line[1]
for chromosome in chromosome_ls:
    print(f'The size of {chromosome} is {chromosome_size_dic[chromosome]}')
    chr_info_dic[chromosome] = 0

a = 0  # a is to moniter the reading status, reaching the end of the file will result in a value of 1
with open(comp_tagged_pair_file, 'r') as pairs_file:
    csv_reader = csv.reader(pairs_file, delimiter='\t')
    line = next(csv_reader)
    with open('/Users/jiangxu/Documents/Seq_data_analysis/dilution_HiC/merge_pair_out/hic_results/data/GM12878_hg19_MQ30/GM12878_diHiC_HindIII_hg19_MQ30_same_comp_contact_ratio.bed', 'w') as bed_file:
        csv_writer = csv.writer(bed_file, delimiter='\t')
        for chromosome in chromosome_ls:   # iterate through the chromosome list to generate bin and collect same compartment contact data
            i = 1 # set the start of the bin to 1
            while ((i+bin_size) < int(chromosome_size_dic[chromosome])) & (a==0) & (chromosome == line[0]): # when bin's upper boarder less than the length of the chromosome and parsing through the file has not reach the end
                if (i+bin_size) <= int(line[1]): # if bin's upper boarder is less than the first pairs read1's position
                    i += bin_size  # then move the bin until the bin cover the read1 of the pair
                    print(f'the last i is{i}')
                    print(f'the last i+bin_size is {i+bin_size}')
                    print(f'the new pair line is {line}')
                    cis_signal_count = 0  # zero cis_count and trans_count when bin move
                    cis_noise_count = 0
                    trans_signal_count = 0
                    trans_noise_count = 0
                    bin_comp = 'N/A'
                for keys in chr_info_dic:
                    chr_info_dic[keys[:]] = 0
                while int(line[1]) <= (i + bin_size):   # when the read1's position is within the bins upper boarder,
                    bin_comp = line[3]
                    r2_comp = line[7]
                    if (line[0] == line[4]) & (bin_comp == r2_comp): # if it is a both cis interaction and a same compartment interaction
                        cis_signal_count += 1
                    if (line[0] != line[4]) & (bin_comp == r2_comp): # if it is a both trans interaction and a same compartment interaction
                        chr_info_dic[line[4]] += 1 # then the same_comp_chr_count will plus 1 to take record of the chromosome
                        trans_signal_count += 1
                    if (line[0] == line[4]) & (bin_comp != r2_comp): # if it is cis and non_same_compartment contact(cis noise)
                        cis_noise_count += 1
                    if (line[0] != line[4]) & (bin_comp != r2_comp): # if it is trans interaction and non_same_compartment contact(trans noise)
                        trans_noise_count += 1
                    try:
                        line = next(csv_reader)
                    except StopIteration:
                        a = 1
                        print(f'a value is {a}')
                        break
                    if chromosome != line[0]:
                        break
                    print(line)

                if (trans_noise_count + trans_signal_count) != 0:
                    trans_same_comp_sig_noise_ratio = round((trans_signal_count / (trans_noise_count + trans_signal_count)), 3)
                if (trans_noise_count + trans_signal_count) == 0:
                    trans_same_comp_sig_noise_ratio = 0
                if (cis_noise_count + cis_signal_count) != 0:
                    cis_same_comp_sig_noise_ratio = round((cis_signal_count /(cis_signal_count + cis_noise_count)), 3)
                if (cis_noise_count + cis_signal_count) == 0:
                    cis_same_comp_sig_noise_ratio = 0
                print(f'the cis signal is {cis_signal_count} for {chromosome} bin {i} to {i+bin_size}')
                print(f'the cis noise is {cis_noise_count} for {chromosome} bin {i} to {i+bin_size}')
                print(f'the trans signal is {trans_signal_count} for {chromosome} bin {i} to {i+bin_size}')
                print(f'the trans noise is {trans_noise_count} for {chromosome} bin {i} to {i+bin_size}')
                print(f'current lower bin number position is {i}')
                chr_info_sorted_ls = sorted(chr_info_dic.items(), key = lambda item: item[1], reverse=True)
                print(f'the sorted same_comp_trans_contact list is {chr_info_sorted_ls}')
                csv_writer.writerow([chromosome, i, i+bin_size, bin_comp, cis_same_comp_sig_noise_ratio, trans_same_comp_sig_noise_ratio, chr_info_sorted_ls[0][0], chr_info_sorted_ls[0][1], chr_info_sorted_ls[1][0], chr_info_sorted_ls[1][1], chr_info_sorted_ls[2][0], chr_info_sorted_ls[2][1]])
            print(f'{chromosome} finished')
#  sorted_same_comp_trans_ls = sorted(chr_info_dic.items(), key=lambda item:item[1], reverse=True)