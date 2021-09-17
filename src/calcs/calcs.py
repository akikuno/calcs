# from parser import parser
from src.calcs import call_cs
import re
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

starts = [int(s.split("\t")[3]) - 1 for s in body]
cigars = [s.split("\t")[5] for s in body]
que_seqs = [s.split("\t")[9] for s in body]

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
ref_seqs = list(compress(ref_fasta, ref_seq_idx)) * len(que_seqs)

###############################################################################
# Trim start sites in reference
###############################################################################


def trim_starts(ref_seq: str, start: int) -> str:
    return ref_seq[start or None:]


ref_seqs_trimed = list(map(trim_starts, ref_seqs, starts))

###############################################################################
# Trim soft/hard clip in query
###############################################################################


def count_clip(cigar):
    _left_clip = re.sub(r'(S|H).*', '', cigar)
    if re.search(r"[A-Z]", _left_clip):
        left_clip_length = 0
    else:
        left_clip_length = int(_left_clip)
    _right_clip = re.sub(r'.*[A-Z]([0-9]+)(S|H)$', r"\1", cigar)
    if re.search(r"[A-Z]", _right_clip):
        right_clip_length = 0
    else:
        right_clip_length = int(_right_clip)
    return left_clip_length, right_clip_length


clip_length = list(map(count_clip, cigars))

que_seqs_clipped = []
append = que_seqs_clipped.append

for idx, (seq, clip) in enumerate(zip(que_seqs, clip_length)):
    left, right = clip
    seq = seq[left or None: -right or None]
    append(seq)

que_seqs_clipped

###############################################################################
# Annotate Insertion in reference
###############################################################################

# re.sub("[0-9]+(S|H)", "", "11115H100M51111S")


def annotate_insertion(ref, cigar):
    start_idx: int = 0
    end_idx: int = 0
    annotates = []
    append = annotates.append
    cigar_trim_clip = re.sub("[0-9]+(S|H)", "", cigar)
    for _cigar in re.split("([0-9]+I)", cigar_trim_clip):
        if "I" in _cigar:
            ins_num = int(_cigar.replace("I", ""))
            append("I" * ins_num)
        else:
            cigar_num = int(re.sub("[MDNSHPX=]", "", _cigar))
            end_idx = start_idx + cigar_num
            append(ref[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


# [annotate_insertion(x, y) for x, y in zip(ref_seqs, cigars)]

ref_seqs_anno = list(map(annotate_insertion, ref_seqs_trimed, cigars))

###############################################################################
# Annotate Deletion in query
###############################################################################


def annotate_deletion(que, cigar):
    start_idx: int = 0
    end_idx: int = 0
    annotates = []
    append = annotates.append
    cigar_trim_clip = re.sub("[0-9]+(S|H)", "", cigar)
    for _cigar in re.split("([0-9]+D)", cigar_trim_clip):
        if "D" in _cigar:
            del_num = int(_cigar.replace("D", ""))
            append("D" * del_num)
        else:
            cigar_num = int(re.sub("[MINPX=]", "", _cigar))
            end_idx = start_idx + cigar_num
            append(que[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


que_seqs_anno = list(map(annotate_deletion, que_seqs_clipped, cigars))

###############################################################################
# Calculate CS tag
###############################################################################

cstags = list(map(call_cs.call_cs_long, ref_seqs_anno, que_seqs_anno))

[print(s) for s in cstags]
if args.long == False:
    cstags = list(map(call_cs.call_cs_short, cstags))


###############################################################################
# Format SAM or PAF
###############################################################################


# if __name__ == "__main__":
#     # query, reference, long, paf, threads = parser()
#     print(args.query)
#     print(args.reference)
#     print(args.long)
#     print(args.paf)
#     print(args.threads)
