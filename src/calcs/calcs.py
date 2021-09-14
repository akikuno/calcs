import pysam
import argparse
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs='?', type=argparse.FileType(),
                        default=sys.stdin)
    args = parser.parse_args()

    file = []
    with args.infile as f:
        for s in f:
            row = s.strip().split("\t")
            file.append(row)
        print(file)


if __name__ == "__main__":
    main()

# parser = argparse.ArgumentParser(
#     description='このプログラムの説明（なくてもよい）')    # 2. パーサを作る

# for x in dir(pysam):
#     print(x, ':', type(eval("pysam."+x)))

# samfile = pysam.AlignmentFile("tests/deletion/test_del.sam", "r")
# samfile

# pysam.calmd("tests/deletion/test_del.sam", "tests/deletion/test_ref.fa")

# pysam.sort("-o", "tmp.bam", "tests/deletion/test_del.sam")

# subprocess.run(
#     ["samtools", "calmd", "tests/deletion/test_del.sam", "tests/deletion/test_ref.fa"])

# pysam.calmd
