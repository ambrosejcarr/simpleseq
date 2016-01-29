#!/usr/local/bin/python3

from sys import argv
import simpleseq


def main(samfile, num_reads=int(1e8)):
    simpleseq.sam.get_RMT_histogram(samfile, num_reads)


if __name__ == "__main__":
    main(*argv[1:])