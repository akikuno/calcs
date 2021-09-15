# import argparse
# import sys


# def parser():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('query', nargs="?", type=argparse.FileType("r"),
#                         default=sys.stdin, help="Input SAM file")
#     parser.add_argument("-r", "--reference", required=True,
#                         help="Reference FASTA file")
#     parser.add_argument("-l", "--long", action="store_true",
#                         help="CS tag as long format")
#     parser.add_argument("-p", "--paf", action="store_true",
#                         help="Output PAF format")
#     parser.add_argument("-@", "--threads", default=1, type=int)
#     args = parser.parse_args()

#     return args.query, args.reference, args.long, args.paf, args.threads
