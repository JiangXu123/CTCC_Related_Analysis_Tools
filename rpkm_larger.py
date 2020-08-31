import csv

file_cryomilling = "/Users/jiangxu/Documents/article/Genome_3D/RPKM_data/CTCC_SDS_HindIII_sole_cryomilling_breaking_hg19_rpkm.tsv"
CTCC_HindIII = "/Users/jiangxu/Documents/article/Genome_3D/RPKM_data/CTCC_SDS_HindIII_sole_cryomilling_breaking_hg19_rpkm.tsv"
with open(CTCC_HindIII, 'r') as file1:
    csv_reader = csv.reader(file1, delimiter='\t')
    for line in csv_reader:
        if float(line[6]) > 2:
            print(line)

