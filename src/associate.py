###############################################################################
# Import modules
###############################################################################

# buidin modules
from time import time as tt  # ?-------------------===
import sys
from concurrent.futures import ProcessPoolExecutor
from itertools import compress

# custom modules
import parse
import load
import trim
import annotate
import call_cstag
import convert

###############################################################################
# Associate reference sequence to each query read
###############################################################################

ARGS_REFERENCE = "../tests/data/associate/reference.mfa"
ARGS_QUERY = open("../tests/data/associate/query.sam")

QUESAM = tuple(load.sam(ARGS_QUERY))

is_header = [_.startswith("@") for _ in QUESAM]
is_alignment = [not _ for _ in is_header]

HEADER = tuple(list(compress(QUESAM, is_header)))
IS_MINIMAP2 = any(["minimap2" in _ for _ in HEADER])
ALIGNMENTS = tuple(list(compress(QUESAM, is_alignment)))

is_mapped = ["*" != s.split("\t")[5] for s in ALIGNMENTS]
is_unmapped = [not _ for _ in is_mapped]

ALIGNMENTS_MAPPED = tuple(list(compress(ALIGNMENTS, is_mapped)))
ALIGNMENTS_UNMAPPED = tuple(list(compress(ALIGNMENTS, is_unmapped)))

RNAMES = tuple([s.split("\t")[2] for s in ALIGNMENTS_MAPPED])
STARTS = tuple([int(s.split("\t")[3]) - 1 for s in ALIGNMENTS_MAPPED])
CIGARS = tuple([s.split("\t")[5] for s in ALIGNMENTS_MAPPED])
QUESEQS = tuple([s.split("\t")[9] for s in ALIGNMENTS_MAPPED])


# Parse Reference FASTA file
dict_fasta = load.fasta(ARGS_REFERENCE)
REFSEQS = tuple([dict_fasta[rname] for rname in RNAMES])
