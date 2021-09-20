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

## Usage

```bash
calcs [options] <in.sam> -r/--reference <in.fasta>
```

## Getting Started

```bash
# CS tag (short form)
calcs aln.sam -r ref.fasta > aln_cs.sam
# CS tag (long form)
calcs aln.sam -r ref.fasta -l > aln_cslong.sam
# PAF format with CS tag (short form)
calcs aln.sam -r ref.fasta -p > aln_cs.paf
# PAF format with CS tag (long form)
calcs aln.sam -r ref.fasta -p -l > aln_cslong.paf
# Multithreading
calcs aln.sam -r ref.fasta -t 4 > aln_cs.sam
```


## Options

```bash
-l, --long: output the cs tag in  the long form 
-t, --threads INT: number of threads to use (default: 1)
```

## Examples

```bash
calcs examples/example.sam -r examples/ref.fa > example_cs.sam
```

If the input file is a BAM/CRAN format, you can use `samtools view`.

```bash
samtools view examples/example.bam |
  calcs -l -r examples/ref.fa |
  samtools sort > example_cslong.bam
```

## `paftools.js sam2paf` vs `calcs`

|                            | sam2paf                    | calcs      |
| -------------------------- | -------------------------- | ---------- |
| Speed                      | +                          | -          |
| Report substitution        | + (if SAM includes MD tag) | +          |
| Report CS tag in long form | + (if SAM includes MD tag) | +          |
| Output format              | PAF                        | SAM or PAF |

## Citation

