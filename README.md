## Description

Append the [minimap2's CS tag](https://github.com/lh3/minimap2#cs) to SAM file.  
If the CS tag is already present, this command will overwrite the existing tag.


> :warning: This tool will be maintained until [the samtools team implements official CS caller](https://github.com/samtools/samtools/issues/1264).

## Installation

You can install `calcs` using pip:

```bash
pip install calcs
```

Alternatively, you can get `calcs` from bioconda:

```
conda install -c bioconda calcs
```

## Getting Started

```bash
calcs aln.sam -r ref.fasta > aln_cs.sam
```


## Options

```bash
-l, --long: Encoding the long form of cs tag (default: false).
-t, --threads INT: Number of threads to use (default: 1).
```

## Examples

```bash
calcs -b -t 2 examples/test.sam -r examples/ref.fa > test_cs.sam
```

If input file is a BAM format, you can use `samtools`.

```bash
samtools view examples/test.bam |
  calcs -l -t 2 -r examples/ref.fa |
  samtools sort > test_cs_long.bam
```

## Citation

