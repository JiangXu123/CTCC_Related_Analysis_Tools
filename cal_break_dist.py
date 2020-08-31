#! /usr/bin/env python

import pandas as pd
import argparse
import time


def run(args):
    start1 = time.perf_counter()
    input_file = args.input
    output_file = args.output
    df = pd.read_csv(input_file, names=['RNAME', 'start', 'breaking_pos', 'count'])
    for index, row in df.iterrows():
        df.at[index, 'dist'] = row[]