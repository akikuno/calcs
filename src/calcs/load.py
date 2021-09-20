import re


def sam(args_query):
    with args_query as f:
        que_sam = [s.strip() for s in f]
    return que_sam


def fasta(args_reference):
    regex = re.compile("(>.*?)\n([ATGCNatgcn\n]*)", re.DOTALL)
    with open(args_reference, 'r') as f:
        content = f.read()
        fasta_dict = {}
        for i in re.findall(regex, content):
            fasta_dict[i[0].replace(">", "")] = i[1].replace('\n', '').upper()
        return fasta_dict
