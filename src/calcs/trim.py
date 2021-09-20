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
