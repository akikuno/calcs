# from parser import parser
from itertools import zip_longest
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

file_que_sam = "tests/subindel/subindel.sam"
file_ref_fasta = "tests/random_100bp.fa"

###############################################################################
# Parse query SAM
###############################################################################

que_sam = []
append = que_sam.append
with open(file_que_sam) as f:
    for s in f:
        row = s.strip()
        append(row)

header_idx = [_.startswith("@") for _ in que_sam]
body_idx = [not _ for _ in header_idx]

header = list(compress(que_sam, header_idx))
body = list(compress(que_sam, body_idx))

start = [int(s.split("\t")[3]) - 1 for s in body]
cigar = [s.split("\t")[5] for s in body]
que_seq = [s.split("\t")[9] for s in body]

###############################################################################
# Parse reference FASTA
###############################################################################

ref_fasta = []
append = ref_fasta.append
with open(file_ref_fasta) as f:
    for s in f:
        row = s.strip()
        append(row)

ref_seq_idx = [not _.startswith(">") for _ in ref_fasta]
ref_seq = list(compress(ref_fasta, ref_seq_idx)) * len(que_seq)

###############################################################################
# Remove soft/hard clip in query sequence
###############################################################################


def count_clip(cigar_str):
    _left_clip = re.sub(r'(S|H).*', '', cigar_str)
    if re.search(r"[A-Z]", _left_clip):
        left_clip_length = 0
    else:
        left_clip_length = int(_left_clip)
    _right_clip = re.sub(r'.*[A-Z]([0-9]+)(S|H)$', r"\1", cigar_str)
    if re.search(r"[A-Z]", _right_clip):
        right_clip_length = 0
    else:
        right_clip_length = int(_right_clip)
    return left_clip_length, right_clip_length


clip_length = list(map(count_clip, cigar))

que_seq_clipped = []
append = que_seq_clipped.append

for idx, (seq, clip) in enumerate(zip(que_seq, clip_length)):
    left, right = clip
    seq = seq[left or None: -right or None]
    append(seq)

que_seq_clipped

###############################################################################
# Annotate Insertion in reference
###############################################################################

tmp_ref_seq: str = "AAATTT"
cigar_str: str = "3M2I2D3M"
expected = "AAAIITTT"
tmp_ref_seq[:3] + int(2) * "I" + tmp_ref_seq[3:6]

# cigar_split = re.split("([0-9]+I)", cigar_str)

ref_seq_ins_annotated = []
append = ref_seq_ins_annotated.append
start_idx: int = 0
for _cigar in cigar_split:
    if "S" in _cigar or "H" in _cigar:
        pass
    elif "I" in _cigar:
        ins_num = int(_cigar.replace("I", ""))
        append("I" * ins_num)
    else:
        print(start_idx)
        cigar_num = int(re.sub("[MDNSHPX=]", "", _cigar))
        end_idx = start_idx + cigar_num
        append(tmp_ref_seq[start_idx or None:end_idx or None])
        start_idx += cigar_num

ref_seq_ins_annotated = ''.join(ref_seq_ins_annotated)

###############################################################################
# Annotate Deletion in query
###############################################################################
cigar
que_seq
re.sub(r"[A-Z]([0-9]+)I", r"\1", cigar_str)

_left_clip = re.sub(r'(S|H).*', '', cigar_str)
if re.search(r"[A-Z]", _left_clip):
    left_clip_length = 0
else:
    left_clip_length = int(_left_clip)
_right_clip = re.sub(r'.*[A-Z]([0-9]+)(S|H)$', r"\1", cigar_str)
if re.search(r"[A-Z]", _right_clip):
    right_clip_length = 0
else:
    right_clip_length = int(_right_clip)
return left_clip_length, right_clip_length


clip_length = list(map(count_clip, cigar))


ref_seq_clipped = []
append = ref_seq_clipped.append

for s, que, ref in zip(start, que_seq_clipped, ref_seq):
    _ = ref[s or None:]
    append(_[:len(que)])

###############################################################################
# Add Insertion in reference
###############################################################################

que_seq_clipped[2]
ref_seq_clipped[2]

# if __name__ == "__main__":
#     # query, reference, long, paf, threads = parser()
#     print(args.query)
#     print(args.reference)
#     print(args.long)
#     print(args.paf)
#     print(args.threads)
