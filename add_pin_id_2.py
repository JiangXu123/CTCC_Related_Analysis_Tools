import argparse
import tabix
import pandas as pd
import concurrent.futures


input_pixel_file = args.input  #'/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/block0.tmp'
bin_f = '/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/hg38_bined_chrom_200b.bed.gz'
output_file = '/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/block0_edited.tmp'

tb = tabix.open(bin_f)  # load tabix with the modified bin file
df = pd.read_csv(input_pixel_file, delimiter='\t')
df.columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'contact_fq']
df['bin1_id'] = None
df['bin2_id'] = None
n = 5000000
split_data = [df[i:i+n] for i in range(0, 20000000, n)]


def add_bins(d):
    for index, row in d.iterrows():    # tabix will return a iterator object, the content can only be obtained via iteration
        results1 = tb.query(row['chrom1'], row['start1'], row['end1'])
        results2 = tb.query(row['chrom2'], row['start2'], row['end2'])
        d.at[index, 'bin1_id'] = int(next(results1)[3])
        d.at[index, 'bin2_id'] = int(next(results2)[3])
    new_data = d.drop(columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
    new_data = new_data[['bin1_id', 'bin2_id', 'contact_fq']]
    return new_data


with concurrent.futures.ProcessPoolExecutor() as executor:
    results = [executor.submit(add_bins, data) for data in split_data]
    for f in concurrent.futures.as_completed(results):
        print(f.result().head())




# new_df = pd.DataFrame(columns=['bin1_id', 'bin2_id', 'contact_fq'])   # create a dataframe with headers to contain new data.
#
# with concurrent.futures.ProcessPoolExecutor() as executor:
#     results = [executor.submit(add_bins, df[i:i+n]) for i in range(0, 20000000, n)]
#     for f in concurrent.futures.as_completed(results):
#         new_df = pd.concat([new_df, f.result()], sort=False)
#
# new_df.to_csv(output_file, index=False, header=True, sep='\t')
