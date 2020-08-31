#! /usr/bin/env python


import pandas as pd
import tabix

in_filename = "/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/CTCC_GM12878_all_394M_hg38_200b_head1M.tsv"
df = pd.read_csv(in_filename, delimiter='\t')
df_trans = df.loc[df['chrom1'] != df['chrom2']]
new_list = pd.DataFrame()
for index, row in df_trans.iterrows():
    if(row[6] > 1):
        new_list = new_list.append(row)
        print(row[6])
print(new_list['count'])
new_list = new_list[['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count']]
new_list.to_csv('/Users/jiangxu/Documents/higlass_data/CTCC_SDS_HindIII/hg38/Trans_over1.tsv', index=False, sep='\t')





