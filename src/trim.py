import re

###############################################################################
# Trim start sites in reference
###############################################################################


def starts(ref_seq: str, start: int) -> str:
    return ref_seq[start or None:]

###############################################################################
# Trim soft-clip in query
###############################################################################


def get_softclip_lengths(cigar):
    """Get the length of left and right softclips"""
    _left = re.sub(r'S.*', '', cigar)
    if re.search(r"[A-Z]", _left):
        left_length = 0
    else:
        left_length = int(_left)
    _right = re.sub(r'.*[A-Z]([0-9]+)S$', r"\1", cigar)
    if re.search(r"[A-Z]", _right):
        right_length = 0
    else:
        right_length = int(_right)
    return left_length, right_length


def softclips(queseq, len_clip):
    left, right = len_clip
    return queseq[left or None: -right or None]
