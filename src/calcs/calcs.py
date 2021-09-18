# from parser import parser
from math import log2
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

args_long = False
args_paf = True
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

###############################################################################
# Annotate Insertion in reference
###############################################################################


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
            cigar_num = sum([int(s or 0)
                            for s in re.split("[MDNPX=]", _cigar)])
            end_idx = start_idx + cigar_num
            append(ref[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


ref_seqs_anno = list(map(annotate_insertion, ref_seqs_trimed, cigars))

###############################################################################
# Annotate Deletion in query
###############################################################################

que = que_seqs_clipped[9]
cigar = cigars[9]


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
            cigar_num = sum([int(s or 0)
                            for s in re.split("[MINPX=]", _cigar)])
            end_idx = start_idx + cigar_num
            append(que[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


que_seqs_anno = list(map(annotate_deletion, que_seqs_clipped, cigars))

###############################################################################
# Calculate CS tag
###############################################################################


def call_cs_long(ref: str, que: str):
    cslong = []
    append = cslong.append
    _cs: str = ''
    _previous: str = ''
    for _ref, _que in zip(list(ref), list(que)):
        # Match
        if _ref == _que and _previous == "M":
            _cs = _ref
        elif _ref == _que and not _previous == "M":
            _cs = "=" + _ref
            _previous = "M"
        # Deletion
        elif _que == "D" and _previous == "D":
            _cs = _ref.lower()
        elif _que == "D" and not _previous == "D":
            _cs = "-" + _ref.lower()
            _previous = "D"
        # Insertion
        elif _ref == "I" and _previous == "I":
            _cs = _que.lower()
        elif _ref == "I" and not _previous == "I":
            _cs = "+" + _que.lower()
            _previous = "I"
        # Substitution
        elif _ref != _que:
            _cs = "*" + _ref.lower() + _que.lower()
            _previous = "S"
        append(_cs)
    return "cs:Z:" + ''.join(cslong)


def call_cs_short(cslong):
    cs = []
    append = cs.append
    cs_split = re.split("(=[A-Z]+)", cslong)
    for _cs in cs_split:
        if _cs.startswith("="):
            _cs = ":" + str(len(re.findall("[A-Z]", _cs)))
        append(_cs)
    return ''.join(cs)


cstags = list(map(call_cs_long, ref_seqs_anno, que_seqs_anno))


if not args_long:
    cstags = list(map(call_cs_short, cstags))


###############################################################################
# Output SAM or PAF
###############################################################################

# PAF format:
# https://github.com/lh3/miniasm/blob/master/PAF.md

if args_paf:
    def determine_strand(flag: int) -> str:
        _cache = flag
        _strand = "+"
        while _cache > 0:
            _power = int(log2(_cache))
            _cache = _cache - 2**_power
            if _power == 4:
                _strand = "-"
                break
        return _strand

    def convert_to_paf(body, cstag, clip_length, start, ref_seq, ref_seq_anno) -> str:
        body_split = body.split("\t")
        start_clips, end_clips = clip_length
        _quename = body_split[0]
        _quelen = str(len(body_split[9]))
        _questart = str(start_clips)
        _queend = str(len(body_split[9]) - end_clips)
        _strand = determine_strand(int(body_split[1]))
        _refname = body_split[2]
        _reflen = str(len(ref_seq))
        _refstart = str(start)
        _refend = str(len(ref_seq_anno.replace("I", "")) + start)
        if "cs:Z:=" in cstag:
            _matches = len(re.findall('[A-Z]', cstag)) - 1
        else:
            cs_split = re.split(r"[a-zA-Z:\-+*]", cstag)
            _matches = sum([int(s or 0) for s in cs_split])
        _matches = str(_matches)
        _blocklen = str(len(ref_seq_anno))
        _quality = body_split[4]
        _others = body_split[11:]
        _others.append(cstag)
        paf = [_quename, _quelen, _questart, _queend, _strand, _refname,
               _reflen, _refstart, _refend, _matches, _blocklen, _quality]
        return '\t'.join(paf + _others)

    paf_cstag = list(map(convert_to_paf, body, cstags, clip_length,
                         starts, ref_seqs, ref_seqs_anno))
    sys.stdout.write('\n'.join(paf_cstag + [""]))

else:
    def insert_cstag(body: str, cstag: str) -> str:
        return re.sub("(\trl:i:0)", "\t" + cstag + r"\1", body)
    body_cstag = list(map(insert_cstag, body, cstags))
    que_sam_cstag = header + body_cstag
    sys.stdout.write('\n'.join(que_sam_cstag + [""]))


# if __name__ == "__main__":
#     # query, reference, long, paf, threads = parser()
#     print(args.query)
#     print(args.reference)
#     print(args.long)
#     print(args.paf)
#     print(args.threads)
