###############################################################################
# Annotate Insertion in reference
###############################################################################
import re

# tmp_ref_seq: str = "AAATTT"
# cigar_str: str = "3M2I2D3M"
# expected = "AAAIITTT"
# tmp_ref_seq[:3] + int(2) * "I" + tmp_ref_seq[3:6]


def annotate_insertion_to_reference(ref, cigar):
    start_idx: int = 0
    end_idx: int = 0
    annotates = []
    append = annotates.append
    for _cigar in re.split("([0-9]+I)", cigar):
        if "S" in _cigar or "H" in _cigar:
            pass
        elif "I" in _cigar:
            ins_num = int(_cigar.replace("I", ""))
            append("I" * ins_num)
        else:
            cigar_num = int(re.sub("[MDNSHPX=]", "", _cigar))
            end_idx = start_idx + cigar_num
            append(ref[start_idx or None:end_idx or None])
            start_idx += cigar_num
    return ''.join(annotates)


ref_seq_ins_annotated = annotate_insertion_to_reference(tmp_ref_seq, cigar_str)
