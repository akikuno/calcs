import re


###############################################################################
# Annotate Insertion in reference
###############################################################################


def insertion(ref, cigar):
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
            cigar_num = sum(
                [int(s or 0) for s in re.split("[MDNPX=]", _cigar)]
            )
            end_idx = start_idx + cigar_num
            append(ref[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


###############################################################################
# Annotate Deletion in query
###############################################################################


def deletion(que, cigar):
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
