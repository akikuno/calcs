# from parser import parser
from itertools import compress
import formatter
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


query = []
with args.query as f:
    for s in f:
        row = s.strip()
        query.append(row)

header = [s.startswith("@") for s in query]
print(list(compress(query, header)))

# print(query)

if __name__ == "__main__":
    # query, reference, long, paf, threads = parser()
    print(args.query)
    print(args.reference)
    print(args.long)
    print(args.paf)
    print(args.threads)
