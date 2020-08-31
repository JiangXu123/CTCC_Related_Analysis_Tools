import cooler
import concurrent.futures
import time


start = time.perf_counter()

in_filename = '/home/rcf-40/jiangxu/proj3/NGS_Sequencing/20180730_CTCC_cryomilling_and_crosslinking/data/all_out_hg38_merged/hic_results/data/CTCC_GM12878_all/pairs_file/mq30_pairs/CTCC_GM12878_all_hg38_MQ30_200b.cool'
read_size = 20000000
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
    processes = [executor.submit(calculate_frequency, cs, (cs + read_size)) for cs in range(0, total_rows, read_size)]
    num = 0
    for f in concurrent.futures.as_completed(processes):
        f.result().to_csv(f'/home/rcf-40/jiangxu/proj3/NGS_Sequencing/20180730_CTCC_cryomilling_and_crosslinking/data/all_out_hg38_merged/hic_results/data/CTCC_GM12878_all/pairs_file/mq30_pairs/block{num}.tmp', header=False, index=False)
        num += 1

end = time.perf_counter()

print(f'multiprocessing uses {end - start} secs')
