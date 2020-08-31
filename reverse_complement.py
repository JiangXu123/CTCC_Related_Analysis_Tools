#! /usr/bin/env python

import argparse


def run(args):
    action_sel = args.action
    input_string = args.input
    compl_bases = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
    if action_sel == 'r':
        ls = list(input_string)
        ls.reverse()
        print(''.join(ls))
    if action_sel == 'rc':
        ls = list(input_string)
        ls.reverse()
        print(''.join([compl_bases[base] for base in ls]))
    if action_sel == 'c':
        print(''.join([compl_bases[base] for base in list(input_string)]))


def main():
    parser = argparse.ArgumentParser(description="reverse or reverse complement a DNA sequence")
    parser.add_argument("-act", help="reverse complement(rc) or reverse(r) or complementary(c)", dest="action", type=str, required=True)
    parser.add_argument("-in", help="input DNA sequence in upper case", dest="input", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()