# from parser import parser
from itertools import compress
# import formatter
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('query', nargs="?", type=argparse.FileType("r"),
                    default=sys.stdin,
                    help="Give the full path to a SAM file")
parser.add_argument("-r", "--reference", required=True,
                    help="Give the full path to a reference FASTA file")
parser.add_argument("-l", "--long", action="store_true",
                    help="Output the cs tag in the long form")
parser.add_argument("-p", "--paf", action="store_true",
                    help="Output PAF")
parser.add_argument("-@", "--threads", default=1, type=int,
                    help="Number of threads [default: 1]")
args = parser.parse_args()

# *TEST========================================================
# arg.query=open("tests/substitution/sub.sam")
query_sam = []
with open("tests/deletion/del.sam") as f:
    for s in f:
        row = s.strip()
        query_sam.append(row)

header = [_.startswith("@") for _ in query_sam]
body = [not _ for _ in header]

header = list(compress(query_sam, header))
body = list(compress(query_sam, body))

cigar = [s.split("\t")[5] for s in body]
query = [s.split("\t")[9] for s in body]


# *TEST========================================================

# print(query)

if __name__ == "__main__":
    # query, reference, long, paf, threads = parser()
    print(args.query)
    print(args.reference)
    print(args.long)
    print(args.paf)
    print(args.threads)
