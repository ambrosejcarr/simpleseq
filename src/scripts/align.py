#!/usr/local/bin/python3

import argparse
import simpleseq
import multiprocessing

def parse_args():
    p = argparse.ArgumentParser(description='merge annotations from drop-seq and in-drop')
    p.add_argument('-m', '--merged-fastq', nargs='+', metavar='G',
                   help='fastq file(s) containing annotated genomic information')
    p.add_argument('-i', '--index', metavar='I', help='genome index folder')
    p.add_argument('-o', '--output-directory', metavar='O', help='name of output directory')
    return p.parse_args()

if __name__ == '__main__':
    args = parse_args()
    n_threads = multiprocessing.cpu_count() - 1
    simpleseq.sam.STAR.align(fastq_file=args.merged_fastq, index=args.index, n_threads=n_threads,
                             output_prefix=args.output_directory)

