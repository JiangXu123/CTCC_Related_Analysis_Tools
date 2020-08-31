import csv

input_file = "/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/ENCFF833FTF_sorted.bed"
output_file = '/Users/jiangxu/Documents/article/Genome_3D/Pileup_assay/ENCFF833FTF_sorted_4.7above.bed'
threshold = 4.7
with open(input_file, 'r') as file1:
    with open(output_file, 'w') as file2:
        csv_reader = csv.reader(file1, delimiter='\t')
        csv_writer = csv.writer(file2, delimiter='\t')
        for line in csv_reader:
            if float(line[8]) > threshold:
                csv_writer.writerow([line[0], line[1], line[2]])
