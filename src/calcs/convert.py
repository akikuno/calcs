import re
from math import log2

# ? PAF format: # ?
# ? https://github.com/lh3/miniasm/blob/master/PAF.md # ?


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


def to_paf(alignment, cstag, len_clip, start, reflen, refseq_anno) -> str:
    alignment_fields = alignment.split("\t")
    start_clips, end_clips = len_clip
    _quename = alignment_fields[0]
    _quelen = str(len(alignment_fields[9]))
    _questart = str(start_clips)
    _queend = str(len(alignment_fields[9]) - end_clips)
    _strand = determine_strand(int(alignment_fields[1]))
    _refname = alignment_fields[2]
    _reflen = str(reflen)
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


def insert_cstag(alignment: str, cstag: str, is_minimap2: bool) -> str:
    if is_minimap2:
        align_with_cstag = alignment.replace("rl:i:", cstag + "\trl:i:")
    else:
        align_with_cstag = alignment + "\t" + cstag
    return align_with_cstag
