#! /usr/bin/env python

import argparse
import tabix
import pandas as pd
import concurrent.futures


def run(args):
    input_pixel_file = args.input
    output_file = args.output
    bin_f = args.bin_file
    n = args.chunk
    tb = tabix.open(bin_f)  # load tabix with the modified bin file
    df = pd.read_csv(input_pixel_file, delimiter='\t')
    df_length = len(df)
    print(df_length)
    def add_bins(data):
        data.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'contact_fq']
        data['bin1_id'] = None
        data['bin2_id'] = None
        for index, row in data.iterrows():    # tabix will return a iterator object, the content can only be obtained via iteration
            results1 = tb.query(row['chrom1'], row['start1'], row['end1'])
            results2 = tb.query(row['chrom2'], row['start2'], row['end2'])
            data.at[index, 'bin1_id'] = int(next(results1)[3])
            data.at[index, 'bin2_id'] = int(next(results2)[3])
        data = data.drop(columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
        return data

    new_df = pd.DataFrame(columns=['bin1_id', 'bin2_id', 'contact_fq'])   # create a dataframe with headers to contain new data.
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(add_bins, df[i:i+n]) for i in range(0, 20000000, 2000000)]
        for f in concurrent.futures.as_completed(results):
            new_df = pd.concat([new_df, f.result()], sort=False)
    new_df.to_csv(output_file, index=False, header=True, sep='\t')


def main():
    parser = argparse.ArgumentParser(description="transform cooler dump join pixel file into simple form with bin_ID")
    parser.add_argument("-in", help="input pairs file", dest="input", type=str, required=True)
    parser.add_argument("-out", help="output files name", dest="output", type=str, required=True)
    parser.add_argument("-mbin", help="modified bin compressed bin file", dest="bin_file", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()