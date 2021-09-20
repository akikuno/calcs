def sam(args_query):
    que_sam = []
    append = que_sam.append
    with args_query as f:
        for s in f:
            row = s.strip()
            append(row)
    return que_sam


def fasta(args_reference):
    ref_fasta = []
    append = ref_fasta.append
    with open(args_reference) as f:
        for s in f:
            row = s.strip()
            append(row)
    return ref_fasta
