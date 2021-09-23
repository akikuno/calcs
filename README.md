## Description

Append the [minimap2's CS tag](https://github.com/lh3/minimap2#cs) to SAM file.  


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

## Usage

```bash
calcs [options] <in.sam> -r/--reference <in.fasta>
```

## Options

```bash
-l/--long: output the cs tag in the long form
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

[paftools.js sam2paf](https://github.com/lh3/minimap2/blob/master/misc/README.md) is anoter command to generate CS tag.  
Here is the comparison between `sam2paf` and `calcs`.

|                     | sam2paf                    | calcs      |
| ------------------- | -------------------------- | ---------- |
| Speed               | +                          | -          |
| Report substitution | + (if SAM includes MD tag) | +          |
| Report CS tag       | + (if SAM includes MD tag) | +          |
| Output format       | PAF                        | SAM or PAF |


