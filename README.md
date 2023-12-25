[![licence](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://choosealicense.com/licenses/mit/)
[![PyPI version](https://img.shields.io/badge/Install%20with-PyPI-brightgreen.svg?style=flat-square)](https://pypi.org/project/calcs/)
[![install with bioconda](https://img.shields.io/badge/Install%20with-Bioconda-brightgreen.svg?style=flat-square)](https://anaconda.org/bioconda/calcs)

## Description

`calcs` is a command-line tool specifically designed to append a [minimap2's cs tag](https://github.com/lh3/minimap2#cs) to a SAM file.  

> [!CAUTION]  
> # `calcs` is deprecated.
> If your SAM file has MD tags, we recommend using [`cstag-cli`](https://github.com/akikuno/cstag-cli) or [`paftools.js sam2paf`](https://github.com/lh3/minimap2/blob/master/misc/README.md).  
> Even if your SAM/BAM files do not have MD tags, **we recommend using [`samtools calmd`](https://www.htslib.org/doc/samtools-calmd.html) to add MD tags, and `cstag-cli` to add cs tags** since `calcs` requires computational time.  
> See [the *Comparison with other tools* section](https://github.com/akikuno/calcs?tab=readme-ov-file#comparison-with-other-tools) for details.


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
# cs tag (short form)
calcs examples/example.sam --reference examples/ref.fa > example_cs.sam

# cs tag (long form)
calcs examples/example.sam --reference examples/ref.fa --long > example_cslong.sam

# PAF format with cs tag (short form)
calcs examples/example.sam --reference examples/ref.fa --paf > example_cs.paf

# PAF format with cs tag (long form)
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

## Comparison with other tools

Here is the brief comparison between `calcs`, `cstag-cli`, and `sam2paf`.

|               | `calcs`         |  `cstag-cli`                |  `sam2paf`            |
| ------------- | ------------- | ------------------------- | ------------------- |
| Input         | SAM and FASTA |  SAM or BAM with a MD tag |  SAM with a MD tag  |
| Output format |  SAM or PAF   |  SAM or BAM               |  PAF                |
| Speed         | Slow🐢             |  Fast🐇                        |  Fast🐇                  |


Instead of `calcs`, we recommend using either `cstag-cli` or `sam2paf`.  
If a SAM file lacks MD tags, you can first add these tags using `samtools calmd` and then apply `cstag-cli`.  
