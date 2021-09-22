import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('query', nargs="?",
                        type=argparse.FileType("r"),
                        default=sys.stdin,
                        metavar="<in.sam>",
                        help="Give the full path to a SAM file")
    parser.add_argument("-r", "--reference", required=True,
                        metavar="<in.fasta>",
                        help="Give the full path to a reference FASTA file")
    parser.add_argument("-l", "--long",
                        action="store_true",
                        help="Output the cs tag in the long form")
    parser.add_argument("-p", "--paf",
                        action="store_true",
                        help="Output PAF")
    parser.add_argument("-t", "--threads",
                        default=1,
                        type=int,
                        metavar="threads",
                        help="Number of threads [default: 1]",)
    parser.add_argument('-v', '--version', action='version', version='0.0.1')
    args = parser.parse_args()
    os_cpus = int(os.cpu_count())  # len(os.sched_getaffinity(0))
    if args.threads > os_cpus:
        threads = os_cpus
    else:
        threads = args.threads
    return args.query, args.reference, args.long, args.paf, threads
