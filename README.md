[![licence](https://anaconda.org/bioconda/calcs/badges/license.svg)](https://opensource.org/search/node/MIT)
[![conda](https://anaconda.org/bioconda/calcs/badges/installer/conda.svg)](https://anaconda.org/bioconda/calcs)
[![PyPI version](https://badge.fury.io/py/calcs.svg)](https://pypi.org/project/calcs/)

## Description

`calcs` is a command-line tool to append a [minimap2's CS tag](https://github.com/lh3/minimap2#cs) to a SAM file.  


> :warning: This tool will be maintained until [the samtools team implements the official CS tag caller](https://github.com/samtools/samtools/issues/1264).

## Installation

You can install `calcs` using pip:

```bash
pip install calcs
```

Alternatively, you can get `calcs` from bioconda:

```
conda install -c bioconda calcs
```

## Usage

```bash
calcs [options] <in.sam> -r/--reference <in.fasta>
```

## Options

```bash
-l/--long: output a cs tag in the long form
-t/--threads [INT]: number of threads to use (default: 1)
```

## Examples

```bash
# CS tag (short form)
calcs examples/example.sam --reference examples/ref.fa > example_cs.sam

# CS tag (long form)
calcs examples/example.sam --reference examples/ref.fa --long > example_cslong.sam

# PAF format with CS tag (short form)
calcs examples/example.sam --reference examples/ref.fa --paf > example_cs.paf

# PAF format with CS tag (long form)
calcs examples/example.sam --reference examples/ref.fa --paf --long > example_cslong.paf

# Multiprocessing
calcs examples/example.sam --reference examples/ref.fa --threads 4 > example_cs.sam
```

If the input file is a BAM/CRAN format, you can use `samtools view`.

```bash
samtools view examples/example.bam |
  calcs -l -r examples/ref.fa |
  samtools sort > example_cslong.bam
```

## `paftools.js sam2paf` vs `calcs`

[paftools.js sam2paf](https://github.com/lh3/minimap2/blob/master/misc/README.md) is anoter command to generate a CS tag.  
Here is the comparison between `sam2paf` and `calcs`.

|                     | sam2paf                    | calcs      |
| ------------------- | -------------------------- | ---------- |
| Speed               | +                          | -          |
| Call a CS tag       | + (if SAM includes a MD tag) | +          |
| Output format       | PAF                        | SAM or PAF |


