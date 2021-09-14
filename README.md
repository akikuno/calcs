# calcs - Generate the CS tag from SAM/BAM/CRAM

## SYNOPSIS

```bash
calcs [-Eeubr] [-C capQcoef] aln.bam ref.fasta
```

## DESCRIPTION
Generate the CS tag. If the CS tag is already present, this command will give a warning if the CS tag generated is different from the existing tag. Output SAM by default.


## OPTIONS

```bash
-b: Output compressed BAM

-t, --threads INT: Number of input/output compression threads to use in addition to main thread [1].
```

## EXAMPLES
Dump BAQ applied alignment for other SNP callers:

```bash
calcs -b -t 2 aln.bam > aln.cs.bam
```
