###############################################################################
# Import modules
###############################################################################

import argparse
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import compress
from math import log2

###############################################################################
# Argument parse
###############################################################################


def parse_args():
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
    if args.threads > len(os.sched_getaffinity(0)):
        threads = len(os.sched_getaffinity(0))
    else:
        threads = args.threads
    return args.query, args.reference, args.long, args.paf, threads


###############################################################################
# Parse query SAM
###############################################################################


def read_query_sam(args_query):
    que_sam = []
    append = que_sam.append
    with args_query as f:
        for s in f:
            row = s.strip()
            append(row)
    return que_sam

###############################################################################
# Parse reference FASTA
###############################################################################


def read_reference_fasta(args_reference):
    ref_fasta = []
    append = ref_fasta.append
    with open(args_reference) as f:
        for s in f:
            row = s.strip()
            append(row)
    return ref_fasta


###############################################################################
# Trim start sites in reference
###############################################################################


def trim_starts(ref_seq: str, start: int) -> str:
    return ref_seq[start or None:]

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


def trim_queseqs(queseq, len_clip):
    left, right = len_clip
    return queseq[left or None: -right or None]

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

###############################################################################
# Output SAM or PAF
###############################################################################
# PAF format:
# https://github.com/lh3/miniasm/blob/master/PAF.md


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


def insert_cstag(alignment: str, cstag: str) -> str:
    return re.sub("(\trl:i:0)", "\t" + cstag + r"\1", alignment)


###############################################################################
# MAIN
###############################################################################


if __name__ == "__main__":
    # Argument parse
    ARGS_QUERY, ARGS_REFERENCE, ARGS_LONG, ARGS_PAF, ARGS_THREADS = parse_args()

    # Parse Query SAM file
    QUESAM = tuple(read_query_sam(ARGS_QUERY))

    is_header = [_.startswith("@") for _ in QUESAM]
    is_alignment = [not _ for _ in is_header]

    HEADER = tuple(list(compress(QUESAM, is_header)))
    ALIGNMENTS = tuple(list(compress(QUESAM, is_alignment)))

    STARTS = tuple([int(s.split("\t")[3]) - 1 for s in ALIGNMENTS])
    CIGARS = tuple([s.split("\t")[5] for s in ALIGNMENTS])
    QUESEQS = tuple([s.split("\t")[9] for s in ALIGNMENTS])

    # Parse Reference FASTA file
    REFFASTA = tuple(read_reference_fasta(ARGS_REFERENCE))
    is_alignment = [not _.startswith(">") for _ in REFFASTA]

    REFSEQS = tuple(list(compress(REFFASTA, is_alignment)) * len(QUESEQS))

    # Trim start sites in reference sequence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(trim_starts, REFSEQS, STARTS))
        REFSEQS_TRIMMED = tuple(_)

    # Trim soft/hard clip into query sequence
    LEN_CLIPS = list(map(len_clips, CIGARS))

    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(trim_queseqs, QUESEQS, LEN_CLIPS))
        QUESEQS_TRIMMED = tuple(_)

    # Annotate Insertion in reference sequence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(annotate_insertion, REFSEQS_TRIMMED, CIGARS))
        REFSEQS_ANNO = tuple(_)

    # Annotate Deletion into query seuqence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(annotate_deletion, QUESEQS_TRIMMED, CIGARS))
        QUESEQS_ANNO = tuple(_)

    # Calculate CS tags
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        cstags = list(executor.map(call_cs_long, REFSEQS_ANNO, QUESEQS_ANNO))

    if ARGS_LONG:
        CSTAGS = tuple(cstags)
    else:
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            _ = list(executor.map(call_cs_short, cstags))
            CSTAGS = tuple(_)

    # Output PAF or SAM
    if ARGS_PAF:
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            paf_cstags = list(executor.map(convert_to_paf,
                                           ALIGNMENTS, CSTAGS, LEN_CLIPS,
                                           STARTS, REFSEQS, REFSEQS_ANNO))
            sys.stdout.write('\n'.join(paf_cstags + [""]))
    else:
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            alignment_cstags = list(executor.map(
                insert_cstag, ALIGNMENTS, CSTAGS))
        SAM_CSTAGS = HEADER + tuple(alignment_cstags)
        sys.stdout.write('\n'.join(SAM_CSTAGS + ("",)))
