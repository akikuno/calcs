import os


# query = []
# with args.infile as f:
#     for s in f:
#         row = s.strip().split("\t")
#         query.append(row)

# reference = []
# with open(args.reference) as f:
#     for s in f:
#         row = s.strip().split("\t")
#         reference.append(row)

# if args.threads > len(os.sched_getaffinity(0)):
#     threads = len(os.sched_getaffinity(0))
# elif args.threads < 1:
#     threads = 1
# else:
#     threads = args.threads
