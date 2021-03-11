import csv

input_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_oqr_second/compartment_pileup/GM12878_subcomp_sorted_Rao_2014.bed"
A1_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/A1.tsv"
A2_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/A2.tsv"
B1_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/B1.tsv"
B2_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/B2.tsv"
B3_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/B3.tsv"
B4_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/B4.tsv"
NA_file = "/Users/jiangxu/Documents/article/Genome_3D/xlink_first_or_second/compartment_pileup/NA.tsv"
with open(input_file, "r") as in_file:
    with open(A1_file, "w") as A1_out:
        with open(A2_file, "w") as A2_out:
            with open(B1_file, "w") as B1_out:
                with open(B2_file, "w") as B2_out:
                    with open(B3_file, "w") as B3_out:
                        with open(B4_file, "w") as B4_out:
                            with open(NA_file, "w") as NA_out:
                                csv_reader = csv.reader(in_file, delimiter="\t")
                                A1_csv_writer = csv.writer(A1_out, delimiter="\t")
                                A2_csv_writer = csv.writer(A2_out, delimiter="\t")
                                B1_csv_writer = csv.writer(B1_out, delimiter="\t")
                                B2_csv_writer = csv.writer(B2_out, delimiter="\t")
                                B3_csv_writer = csv.writer(B3_out, delimiter="\t")
                                B4_csv_writer = csv.writer(B4_out, delimiter="\t")
                                NA_csv_writer = csv.writer(NA_out, delimiter="\t")
                                for line in csv_reader:
                                    if line[3] == "A1":
                                        A1_csv_writer.writerow([line[0], line[1], line[2]])
                                    elif line[3] == "A2":
                                        A2_csv_writer.writerow([line[0], line[1], line[2]])
                                    elif line[3] == "B1":
                                        B1_csv_writer.writerow([line[0], line[1], line[2]])
                                    elif line[3] == "B2":
                                        B2_csv_writer.writerow([line[0], line[1], line[2]])
                                    elif line[3] == "B3":
                                        B3_csv_writer.writerow([line[0], line[1], line[2]])
                                    elif line[3] == "B4":
                                        B4_csv_writer.writerow([line[0], line[1], line[2]])
                                    elif line[3] == "NA":
                                        NA_csv_writer.writerow([line[0], line[1], line[2]])