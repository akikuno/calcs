# calcs - Generate the CS tag from SAM/BAM/CRAM

## Installation

You can install calcs using pip:

```bash
pip install calcs
```

Alternatively, you can get calcs from bioconda:

```
conda install -c bioconda calcs
```

## Getting Started

```bash
calcs aln.sam ref.fasta > aln_cs.sam
```

## Description

Generate the CS tag. If the CS tag is already present, this command will give a warning if the CS tag generated is different from the existing tag.


## Options

```bash
-t, --threads INT: Number of threads to use (default: 1).
```

## Examples

```bash
calcs -b -t 2 examples/test.sam examples/ref.fa > test_cs.sam
```

```bash
# Input BAM file via stdin
samtools view examples/test.bam |
  calcs -b -t 2 examples/ref.fa |
  samtools sort > test_cs.bam
```

## Citation

