#! /usr/bin/env python
# this program is to test the multiprocess function to process_files
import os
import multiprocessing
import argparse


def run(args):
    thread_number = args.t_number
    input_file = args.in_file
    try:
        path = os.path.abspath('.') + '/temp'
        os.makedirs(path)  # make a temporary folder in the current directory
    except:
        print('folder already exist')
    print(path)
    in_ls = []
    out_ls = []
    total_line_count = 0

    for i in range(0, thread_number):
        temp_input_file = './temp/temp_in_' + str(i)
        temp_output_file = './temp/temp_out_' + str(i)
        in_ls.append(temp_input_file)
        out_ls.append(temp_output_file)

    with open(input_file, 'r') as file1:
        for _ in file1:
            total_line_count += 1
    print(total_line_count)
    block_lines = int(total_line_count / thread_number)
    block_index = 0
    block_line_count = 0
    with open(input_file, 'r') as file1:
        for line in file1:
            if block_line_count <= block_lines:
                with open(in_ls[block_index], 'a') as file2:
                    file2.write(line)
            elif block_line_count > block_lines:
                block_line_count = 0
                block_index += 1
                print(f'{block_index} block split finished')
                with open(in_ls[block_index], 'w') as file2:
                    file2.write(line)
            block_line_count += 1

    print(f'file is split to {block_index} blocks into the temp folder')
    # the above code create the temporary folder and splits the file into smaller files
    processes = []

    temp_in_file_name_ls = os.listdir(os.path.abspath('.') + '/temp')
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
    for i in range(0, thread_number):  # the zip method allow to open two list in parallel
        in_file = temp_in_file_path_ls[i]
        out_file = temp_out_file_path_ls[i]
        print(in_file)
        print(out_file)
        proc = multiprocessing.Process(target=process_file, args=[in_file, out_file])
        proc.start()
        processes.append(proc)

    for process in processes:
        process.join()


def process_file(input_file, output_file):
    with open(input_file, 'r') as file1:
        with open(output_file, 'w') as file2:
            for line in file1:
                file2.write(line)


def main():
    parser = argparse.ArgumentParser(description="tagging HiC-Pro pair's sub-compartment")
    parser.add_argument("-in", help="input pairs file", dest="in_file", type=str, required=True)
    parser.add_argument("-tn", help="thread number", dest="t_number", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
