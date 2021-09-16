# from parser import parser
import re
import subprocess
from itertools import compress
# import formatter
import argparse
import sys

# parser = argparse.ArgumentParser()
# parser.add_argument('query', nargs="?", type=argparse.FileType("r"),
#                     default=sys.stdin,
#                     help="Give the full path to a SAM file")
# parser.add_argument("-r", "--reference", required=True,
#                     help="Give the full path to a reference FASTA file")
# parser.add_argument("-l", "--long", action="store_true",
#                     help="Output the cs tag in the long form")
# parser.add_argument("-p", "--paf", action="store_true",
#                     help="Output PAF")
# parser.add_argument("-@", "--threads", default=1, type=int,
#                     help="Number of threads [default: 1]")
# args = parser.parse_args()

# *TEST========================================================
# arg.query=open("tests/deletion/del.sam")
# arg.reference="tests/random_100bp.fa"
# ! posicon
subprocess.run(["cat", "tests/deletion/del.sam"])
subprocess.run(["cat", "tests/deletion/del_cslong.sam"])
###

que_sam = []
with open("tests/deletion/del.sam") as f:
    for s in f:
        row = s.strip()
        que_sam.append(row)

header_idx = [_.startswith("@") for _ in que_sam]
body_idx = [not _ for _ in header_idx]

header = list(compress(que_sam, header_idx))
body = list(compress(que_sam, body_idx))

cigar = [s.split("\t")[5] for s in body]
que_seq = [s.split("\t")[9] for s in body]

ref_fasta = []
with open("tests/random_100bp.fa") as f:
    for s in f:
        row = s.strip()
        ref_fasta.append(row)

ref_seq_idx = [not _.startswith(">") for _ in ref_fasta]
ref_seq = list(compress(ref_fasta, ref_seq_idx))

cigar
que_seq
ref_seq


def cigarsplit(s):
    group_by = 2
    s_split = re.split("(\d+)", s)[1:]
    s_group_by = [s_split[i:i + group_by]
                  for i in range(0, len(s_split), group_by)]
    return s_group_by


[cigarsplit(s) for s in cigar]

# *TEST========================================================

# print(query)

if __name__ == "__main__":
    # query, reference, long, paf, threads = parser()
    print(args.query)
    print(args.reference)
    print(args.long)
    print(args.paf)
    print(args.threads)
