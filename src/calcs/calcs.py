import subprocess
import pysam


for x in dir(pysam):
    print(x, ':', type(eval("pysam."+x)))

samfile = pysam.AlignmentFile("tests/deletion/test_del.sam", "r")
samfile

pysam.calmd("tests/deletion/test_del.sam", "tests/deletion/test_ref.fa")

pysam.sort("-o", "tmp.bam", "tests/deletion/test_del.sam")

subprocess.run(
    ["samtools", "calmd", "tests/deletion/test_del.sam", "tests/deletion/test_ref.fa"])

pysam.calmd
