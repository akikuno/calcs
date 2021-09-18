from math import log2
import re
from itertools import compress
import argparse
import sys

###############################################################################
# Argument parse
###############################################################################

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

args_query = args.query
args_reference = args.reference
args_query = args.query
args_long = args.long
args_paf = args.paf
args_threads = args.threads

# args_query = open("tests/subindel/subindel.sam")
# args_reference = "tests/random_100bp.fa"

###############################################################################
# Parse query SAM
###############################################################################

que_sam = []
append = que_sam.append
with args_query as f:
    for s in f:
        row = s.strip()
        append(row)

is_header = [_.startswith("@") for _ in que_sam]
is_alignment = [not _ for _ in is_header]

HEADER = tuple(list(compress(que_sam, is_header)))
ALIGNMENTS = tuple(list(compress(que_sam, is_alignment)))

STARTS = tuple([int(s.split("\t")[3]) - 1 for s in ALIGNMENTS])
CIGARS = tuple([s.split("\t")[5] for s in ALIGNMENTS])
QUESEQS = tuple([s.split("\t")[9] for s in ALIGNMENTS])

###############################################################################
# Parse reference FASTA
###############################################################################

ref_fasta = []
append = ref_fasta.append
with open(args_reference) as f:
    for s in f:
        row = s.strip()
        append(row)

is_alignment = [not _.startswith(">") for _ in ref_fasta]
REFSEQS = tuple(list(compress(ref_fasta, is_alignment)) * len(QUESEQS))

###############################################################################
# Trim start sites in reference
###############################################################################


def trim_starts(ref_seq: str, start: int) -> str:
    return ref_seq[start or None:]


REFSEQS_TRIMMED = tuple(list(map(trim_starts, REFSEQS, STARTS)))

###############################################################################
# Trim soft/hard clip in query
###############################################################################


def len_clips(cigar):
    """Get the length of left and right clips"""
    _left = re.sub(r'(S|H).*', '', cigar)
    if re.search(r"[A-Z]", _left):
        left_length = 0
    else:
        left_length = int(_left)
    _right = re.sub(r'.*[A-Z]([0-9]+)(S|H)$', r"\1", cigar)
    if re.search(r"[A-Z]", _right):
        right_length = 0
    else:
        right_length = int(_right)
    return left_length, right_length


LEN_CLIPS = list(map(len_clips, CIGARS))

queseqs_trimmed = []
append = queseqs_trimmed.append

for idx, (seq, clip) in enumerate(zip(QUESEQS, LEN_CLIPS)):
    left, right = clip
    seq = seq[left or None: -right or None]
    append(seq)

QUESEQS_TRIMMED = tuple(queseqs_trimmed)

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


refseqs_anno = list(map(annotate_insertion, REFSEQS_TRIMMED, CIGARS))
REFSEQS_ANNO = tuple(refseqs_anno)

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
            cigar_num = sum([int(s or 0)
                            for s in re.split("[MINPX=]", _cigar)])
            end_idx = start_idx + cigar_num
            append(que[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


queseqs_anno = list(map(annotate_deletion, QUESEQS_TRIMMED, CIGARS))
QUESEQS_ANNO = tuple(queseqs_anno)

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


cstags = list(map(call_cs_long, REFSEQS_ANNO, QUESEQS_ANNO))
CSTAGS = tuple(cstags)

if not args_long:
    CSTAGS = tuple(list(map(call_cs_short, cstags)))


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

    def convert_to_paf(alignment, cstag, len_clip, start, refseq, refseq_anno) -> str:
        alignment_fields = alignment.split("\t")
        start_clips, end_clips = len_clip
        _quename = alignment_fields[0]
        _quelen = str(len(alignment_fields[9]))
        _questart = str(start_clips)
        _queend = str(len(alignment_fields[9]) - end_clips)
        _strand = determine_strand(int(alignment_fields[1]))
        _refname = alignment_fields[2]
        _reflen = str(len(refseq))
        _refstart = str(start)
        _refend = str(len(refseq_anno.replace("I", "")) + start)
        if "cs:Z:=" in cstag:
            _matches = len(re.findall('[A-Z]', cstag)) - 1
        else:
            cs_split = re.split(r"[a-zA-Z:\-+*]", cstag)
            _matches = sum([int(s or 0) for s in cs_split])
        _matches = str(_matches)
        _blocklen = str(len(refseq_anno))
        _quality = alignment_fields[4]
        _others = alignment_fields[11:]
        _others.append(cstag)
        paf = [_quename, _quelen, _questart, _queend, _strand, _refname,
               _reflen, _refstart, _refend, _matches, _blocklen, _quality]
        return '\t'.join(paf + _others)

    paf_cstag = list(map(convert_to_paf, ALIGNMENTS, CSTAGS, LEN_CLIPS,
                         STARTS, REFSEQS, REFSEQS_ANNO))
    sys.stdout.write('\n'.join(paf_cstag + [""]))

else:
    def insert_cstag(alignment: str, cstag: str) -> str:
        return re.sub("(\trl:i:0)", "\t" + cstag + r"\1", alignment)
    body_cstag = list(map(insert_cstag, ALIGNMENTS, CSTAGS))
    SAM_CSTAGS = HEADER + tuple(body_cstag)
    sys.stdout.write('\n'.join(SAM_CSTAGS + ("",)))


# if __name__ == "__main__":
#     # query, reference, long, paf, threads = parser()
#     print(args.query)
#     print(args.reference)
#     print(args.long)
#     print(args.paf)
#     print(args.threads)
