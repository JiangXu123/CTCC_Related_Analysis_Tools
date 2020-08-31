#! /usr/bin/env python

import cooler
import argparse
import concurrent.futures
import time


def run(args):
    start = time.perf_counter()

    in_filename = args.input
    # out_filename = args.output
    read_size = args.chsz

    c = cooler.Cooler(in_filename)
    total_contact = c.pixels()[:]['count'].sum()
    total_rows = c.pixels()[:].shape[0]
    print(total_contact)
    print(total_rows)

    def calculate_frequency(s, e):
        df = c.pixels(join=True)[s:e]
        df['count'] = 1000000 * df['count'] / total_contact
        return df

    with concurrent.futures.ProcessPoolExecutor() as executor:
        processes = [executor.submit(calculate_frequency, cs, (cs + read_size)) for cs in
                     range(0, total_rows, read_size)]
        num = 0
        for f in concurrent.futures.as_completed(processes):
            f.result().to_csv(f'/home/rcf-40/jiangxu/proj3/NGS_Sequencing/20180730_CTCC_cryomilling_and_crosslinking/data/all_out_hg38_merged/hic_results/data/CTCC_GM12878_all/pairs_file/mq30_pairs/block{num}.tmp', header=False, index=False)
            num += 1
    end = time.perf_counter()

    print(f'multiprocessing uses {end - start} secs')


def main():
    parser = argparse.ArgumentParser(description="filter_for_high_trans_interaction_bins")
    parser.add_argument("-in", help="input cool file", dest="input", type=str, required=True)
    # parser.add_argument("-out", help="output tsv file", dest="output", type=str, required=False)
    parser.add_argument("-cs", help="chunk size", dest="chsz", type=int, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()