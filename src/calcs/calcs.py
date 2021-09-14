import argparse
import sys
import numpy as np
import pandas as pd
import os


def _argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType("r"),
                        default=sys.stdin)
    parser.add_argument("-r", "--reference", required=True)
    parser.add_argument("-l", "--long", action="store_true")
    parser.add_argument("-t", "--threads", default=1, type=int)
    args = parser.parse_args()

    query = []
    with args.infile as f:
        for s in f:
            row = s.strip().split("\t")
            query.append(row)

    reference = []
    with open(args.reference) as f:
        for s in f:
            row = s.strip().split("\t")
            reference.append(row)

    if args.threads > len(os.sched_getaffinity(0)):
        threads = len(os.sched_getaffinity(0))
    elif args.threads < 1:
        threads = 1
    else:
        threads = args.threads

    return query, reference, args.long, threads


if __name__ == "__main__":
    query, reference, long, threads = _argparse()
    print(query)
    print(reference)
    print(long)
    print(threads)
