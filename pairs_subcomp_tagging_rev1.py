#! /usr/bin/env python
# this is the revision of the program 'pairs_subcomp_tagging.py' to enable multi-processing
import os
import multiprocessing
import argparse
import tabix
import csv
import time
import shutil
import subprocess


def run(args):
    thread_number = args.t_number
    input_file = args.in_file
    output_file = args.out_file
    compartment_filename = args.comp
    experiment = args.exp
    valid_line_count = 0
    total_line_count = 0
    start = time.perf_counter()
    temp_path = os.path.abspath('.') + '/temp'
    try:
        os.makedirs(temp_path)  # make a temporary folder in the current directory
        print(f'The temporary folder created at {temp_path}')
    except:
        print('temp folder already exist')

    with open(input_file, 'r') as file1:
        for _ in file1:
            total_line_count += 1  # count the total number of pairs

    block_lines = int(total_line_count / thread_number)  # make an estimation how many lines should contain per block
    print(f'total of {total_line_count} lines will be split into {block_lines} per temporary block file')
    input_path = os.path.abspath('.') + '/' + input_file
    print(f'input file path is {input_path}')
    os.chdir(temp_path)
    subprocess.call(f'cat {input_path} | split -d -l {block_lines} - temp_in_', shell=True) # I cheated here using bash shell's split method, which is much faster in splitting files
    os.chdir('..')
    block_num = len(os.listdir('temp'))
    print(f'Input file is split to {block_num} blocks into the temp folder')
    # the above code create the temporary folder and splits the file into smaller files
    processes = []
    temp_in_file_name_ls = os.listdir(os.path.abspath('.') + '/' + 'temp')
    temp_in_file_path_ls = []
    temp_out_file_name_ls = []
    temp_out_file_path_ls = []
    for name in temp_in_file_name_ls:
        temp_in_file_path_ls.append(os.path.abspath('.') + '/temp/' + name)
        temp_out_file_name_ls.append('temp_out_' + name.split('_')[-1])

    for name in temp_out_file_name_ls:
        temp_out_file_path_ls.append(os.path.abspath('.') + '/temp/' + name)
    print(temp_in_file_path_ls)
    print(temp_out_file_path_ls)
    for i in range(0, block_num):  # the zip method allow to open two list in parallel
        in_file = temp_in_file_path_ls[i]
        out_file = temp_out_file_path_ls[i]
        print(in_file)
        print(out_file)
        proc = multiprocessing.Process(target=process_file, args=[in_file, compartment_filename, experiment, out_file])
        proc.start()
        processes.append(proc)

    for process in processes:
        process.join()        # wait until all processes finished

    with open(output_file, 'w') as file1:   # concatenate the temporary output files
        for file in temp_out_file_path_ls:
            with open(file, 'r') as file2:
                for line in file2:
                    file1.write(line)
                    valid_line_count += 1
    shutil.rmtree(temp_path)  # remove the temporary folder after finishing processing
    end = time.perf_counter()
    print(f'temporary folder removed')
    print(f'from {total_line_count} pairs, {valid_line_count} pairs tagged in {round(end-start, 2)} sec')


def process_file(in_file, compartment_file, experiment, out_file):
    with open(in_file, 'r') as file1:
        with open(out_file, 'w') as file2:
            csv_reader = csv.reader(file1, delimiter='\t')
            csv_writer = csv.writer(file2, delimiter='\t')
            tb = tabix.open(compartment_file)
            for line in csv_reader:
                try:
                    result1 = tb.query(line[1], int(line[2]) - 1, int(line[2]))
                    result2 = tb.query(line[4], int(line[5]) - 1, int(line[5]))
                    read1_comp = next(result1)[3]
                    read2_comp = next(result2)[3]
                    interaction_type = read1_comp + '-' + read2_comp
                    if line[1] == line[4]:
                        distance = abs(int(line[5]) - int(line[2]))
                    elif line[1] != line[4]:
                        distance = 'NA'
                    csv_writer.writerow([line[1], line[2], line[3], line[4], line[5], line[6], distance, interaction_type, experiment])
                    # columns are: chr_1_name, chr_1_pos, R1_strand, chr_2_name, chr_2_pos, R2_strand, distance, interaction_type, experiment
                except:
                    pass


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-i", help="input pairs file", dest="in_file", type=str, required=True)
    parser.add_argument("-o", help="output compartment tagged pairs file", dest="out_file", type=str, required=True)
    parser.add_argument("-t", help="thread number", dest="t_number", type=int, required=True)
    parser.add_argument("-e", help="experiment", dest="exp", type=str, required=True)
    parser.add_argument("-c", help="tabix indexed compartment file", dest="comp", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
