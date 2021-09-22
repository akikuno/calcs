import re

###############################################################################
# Trim unmapped_region in reference
###############################################################################

p_clip = re.compile(r"[0-9]+(H|S)")
p = re.compile(r"M|D")


def unmapped_region(refseq: str, start: int, cigar: int) -> str:
    '''Trim unmapped_region in reference'''
    _cigar = p_clip.sub("", cigar, count=2)
    _cigar = re.sub(r"[0-9]*I", "", _cigar)
    end = start + sum(map(int, filter(None, p.split(_cigar))))
    return refseq[start: end]


###############################################################################
# Trim soft-clip in query
###############################################################################


def get_softclip_lengths(cigar):
    """Get the length of left and right softclips"""
    left_clips = cigar.split("S", 1)[0]
    try:
        left_clips_len = int(left_clips)
    except ValueError:
        left_clips_len = 0
    right_clips = re.sub(r'.*[A-Z]([0-9]+)S$', r"\1", cigar)
    try:
        right_clips_len = int(right_clips)
    except ValueError:
        right_clips_len = 0
    return left_clips_len, right_clips_len


def softclips(queseq, len_clip):
    left_clips_len, right_clips_len = len_clip
    return queseq[left_clips_len or None: -right_clips_len or None]
