###############################################################################
# Import modules
###############################################################################

# buidin modules
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import compress, chain

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


def main():
    # Argument parse
    ARGS_QUERY, ARGS_REFERENCE, ARGS_LONG, ARGS_PAF, ARGS_THREADS = parse.parse_args()

    # Parse Query SAM file
    QUESAM = tuple(load.sam(ARGS_QUERY))

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
    CHUNKSIZE = int(len(QUESEQS)/ARGS_THREADS)

    # Trim soft clip into query sequence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        LEN_CLIPS = list(executor.map(
            trim.get_softclip_lengths, CIGARS,
            chunksize=CHUNKSIZE)
        )

    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        QUESEQS_TRIMMED = executor.map(
            trim.softclips, QUESEQS, LEN_CLIPS,
            chunksize=CHUNKSIZE
        )

    # Parse Reference FASTA file
    dict_fasta = load.fasta(ARGS_REFERENCE)
    REFSEQS = (dict_fasta[rname] for rname in RNAMES)
    REFLENS = (len(dict_fasta[rname]) for rname in RNAMES)

    # Trim start sites in reference sequence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        REFSEQS_TRIMMED = executor.map(
            trim.unmapped_region, REFSEQS, STARTS, CIGARS,
            chunksize=CHUNKSIZE
        )

    # Annotate Insertion in reference sequence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        REFSEQS_ANNO = list(
            executor.map(
                annotate.insertion, REFSEQS_TRIMMED, CIGARS,
                chunksize=CHUNKSIZE)
        )

    # Annotate Deletion into query seuqence
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        QUESEQS_ANNO = executor.map(
            annotate.deletion, QUESEQS_TRIMMED, CIGARS,
            chunksize=CHUNKSIZE)

    # Calculate CS tags
    with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
        cstags = executor.map(
            call_cstag.long_form, REFSEQS_ANNO, QUESEQS_ANNO,
            chunksize=CHUNKSIZE)

    if ARGS_LONG:
        CSTAGS = cstags
    else:
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            CSTAGS = executor.map(
                call_cstag.short_form, cstags,
                chunksize=CHUNKSIZE
            )
    # Output PAF or SAM
    if ARGS_PAF:
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            paf_cstags = list(
                executor.map(
                    convert.to_paf,
                    ALIGNMENTS_MAPPED, CSTAGS, LEN_CLIPS,
                    STARTS, REFLENS, REFSEQS_ANNO,
                    chunksize=CHUNKSIZE
                )
            )
        try:
            sys.stdout.write('\n'.join(paf_cstags + [""]))
        except (BrokenPipeError, IOError):
            pass
    else:
        with ProcessPoolExecutor(max_workers=ARGS_THREADS) as executor:
            ALIGNMENT_CSTAGS = executor.map(
                convert.insert_cstag,
                ALIGNMENTS_MAPPED,
                CSTAGS,
                [IS_MINIMAP2] * len(ALIGNMENTS_MAPPED),
                chunksize=CHUNKSIZE)

        SAM_CSTAGS = chain(HEADER, ALIGNMENT_CSTAGS, ALIGNMENTS_UNMAPPED)
        try:
            sys.stdout.write('\n'.join(map(str, SAM_CSTAGS)))
        except (BrokenPipeError, IOError):
            pass


###############################################################################
# Call MAIN
###############################################################################

if __name__ == "__main__":
    main()
