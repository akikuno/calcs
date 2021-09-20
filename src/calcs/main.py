###############################################################################
# Import modules
###############################################################################

# buidin modules
from time import time as tt  # ?-------------------===
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import compress

# custom modules
from calcs import parse
from calcs import load
from calcs import trim
from calcs import annotate
from calcs import call_cstag
from calcs import convert

###############################################################################
# MAIN
###############################################################################

# ARGS_QUERY = open("../tests/data/subindel/subindel.sam") # ?-------------------
# ARGS_REFERENCE = "../tests/data/random_100bp.fa" # ?-------------------


def main():
    # Argument parse
    ARGS_QUERY, ARGS_REFERENCE, ARGS_LONG, ARGS_PAF, ARGS_THREADS = parse.parse_args()

    # Parse Query SAM file
    t = tt()  # ?-------------------
    QUESAM = tuple(load.sam(ARGS_QUERY))
    readsam_time = tt() - t  # ?-------------------
    print(f"Read Query: {readsam_time:.2} sec",
          file=sys.stderr)  # ?-------------------

    t = tt()  # ?-------------------
    is_header = [_.startswith("@") for _ in QUESAM]
    HEADER = tuple(list(compress(QUESAM, is_header)))
    IS_MINIMAP2 = any(["minimap2" in _ for _ in HEADER])

    is_alignment = [not _ for _ in is_header]
    alignments = tuple(list(compress(QUESAM, is_alignment)))
    is_mapped = ["*" != s.split("\t")[5] for s in alignments]
    is_unmapped = [not _ for _ in is_mapped]

    ALIGNMENTS_MAPPED = tuple(list(compress(alignments, is_mapped)))
    ALIGNMENTS_UNMAPPED = tuple(list(compress(alignments, is_unmapped)))

    RNAMES = tuple([s.split("\t")[2] for s in ALIGNMENTS_MAPPED])
    STARTS = tuple([int(s.split("\t")[3]) - 1 for s in ALIGNMENTS_MAPPED])
    CIGARS = tuple([s.split("\t")[5] for s in ALIGNMENTS_MAPPED])
    QUESEQS = tuple([s.split("\t")[9].upper() for s in ALIGNMENTS_MAPPED])
    print(f"Parse sam: {tt() - t:.2} sec",
          file=sys.stderr)  # ?-------------------

    # Parse Reference FASTA file
    dict_fasta = load.fasta(ARGS_REFERENCE)
    REFSEQS = tuple([dict_fasta[rname] for rname in RNAMES])

    # Trim start sites in reference sequence
    t = tt()  # ?-------------------
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(trim.starts, REFSEQS, STARTS,
                              chunksize=int(len(REFSEQS)/ARGS_THREADS)))
        REFSEQS_TRIMMED = tuple(_)
    print(f"Trim start sites: {tt() - t:.2} sec",
          file=sys.stderr)  # ?-------------------

    # Trim soft clip into query sequence
    t = tt()  # ?--------------------------------------
    LEN_CLIPS = list(map(trim.get_softclip_lengths, CIGARS))

    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(trim.softclips, QUESEQS, LEN_CLIPS,
                              chunksize=int(len(QUESEQS)/ARGS_THREADS)))
        QUESEQS_TRIMMED = tuple(_)

    print(f"Trim softclips: {tt() - t:.2} sec",
          file=sys.stderr)  # ?-------------------

    # Annotate Insertion in reference sequence
    t = tt()  # ?-------------------
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(annotate.insertion, REFSEQS_TRIMMED, CIGARS))
        REFSEQS_ANNO = tuple(_)
    print(f"Annotate Insertion: {tt() - t:.2} sec",
          file=sys.stderr)  # ?-------------------

    # Annotate Deletion into query seuqence
    t = tt()  # ?-------------------
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        _ = list(executor.map(annotate.deletion, QUESEQS_TRIMMED, CIGARS,
                              chunksize=int(len(QUESEQS_TRIMMED)/ARGS_THREADS)))
        QUESEQS_ANNO = tuple(_)
    print(f"Annotate Deletion: {tt() - t:.2} sec",
          file=sys.stderr)  # ?-------------------

    # Calculate CS tags
    t = tt()  # ?-------------------
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        cstags = list(executor.map(call_cstag.long_form, REFSEQS_ANNO, QUESEQS_ANNO,
                                   chunksize=int(len(REFSEQS_ANNO)/ARGS_THREADS)))
    print(f"Call CSlong: {tt() - t:.2} sec",
          file=sys.stderr)  # ?-------------------

    if ARGS_LONG:
        CSTAGS = tuple(cstags)
    else:
        t = tt()  # ?-------------------
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            _ = list(executor.map(call_cstag.short_form, cstags,
                                  chunksize=int(len(cstags)/ARGS_THREADS)))
            CSTAGS = tuple(_)
        print(f"Call CSshort: {tt() - t:.2} sec",
              file=sys.stderr)  # ?-------------------

    # Output PAF or SAM
    if ARGS_PAF:
        t = tt()  # ?-------------------
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            paf_cstags = list(executor.map(convert.to_paf,
                                           ALIGNMENTS_MAPPED, CSTAGS, LEN_CLIPS,
                                           STARTS, REFSEQS, REFSEQS_ANNO,
                                           chunksize=int(len(ALIGNMENTS_MAPPED)/ARGS_THREADS)))
        print(f"Convert to PAF: {tt() - t:.2} sec",
              file=sys.stderr)  # ?-------------------
        try:
            sys.stdout.write('\n'.join(paf_cstags + [""]))
        except (BrokenPipeError, IOError):
            pass
    else:
        t = tt()  # ?-------------------
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            alignment_cstags = list(executor.map(
                convert.insert_cstag, ALIGNMENTS_MAPPED, CSTAGS,
                [IS_MINIMAP2] * len(ALIGNMENTS_MAPPED),
                chunksize=int(len(ALIGNMENTS_MAPPED)/ARGS_THREADS)))
        print(f"Insert CStag: {tt() - t:.2} sec",
              file=sys.stderr)  # ?-------------------
        SAM_CSTAGS = HEADER + tuple(alignment_cstags) + \
            tuple(ALIGNMENTS_UNMAPPED)
        try:
            sys.stdout.write('\n'.join(SAM_CSTAGS + ("",)))
        except (BrokenPipeError, IOError):
            pass


###############################################################################
# Call MAIN
###############################################################################

if __name__ == "__main__":
    main()
