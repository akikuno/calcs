# calcs - Generate the CS tag from SAM/BAM/CRAM

## SYNOPSIS
calcs [-Eeubr] [-C capQcoef] aln.bam ref.fasta

## DESCRIPTION
Generate the CS tag. If the CS tag is already present, this command will give a warning if the CS tag generated is different from the existing tag. Output SAM by default.

Calmd can also read and write CRAM files although in most cases it is pointless as CRAM recalculates MD and NM tags on the fly. The one exception to this case is where both input and output CRAM files have been / are being created with the no_ref option.

Note that some aligners do not include sequence or confidence values in secondary and supplementary alignment records. Where this happens in SAM files, a “*” character will be seen in the SEQ and QUAL columns. These records will be skipped, as it is not possible to recalculate the MD and NM tags without access to the query sequence. samtools calmd will emit a warning if any records have been skipped for this reason.

## OPTIONS

```
-b
Output compressed BAM

-t, --threads INT
Number of input/output compression threads to use in addition to main thread [1].
```

## EXAMPLES
Dump BAQ applied alignment for other SNP callers:

```bash
calcs -b -t 2 aln.bam > aln.cs.bam
```

## AUTHOR

Written by Akihiro Kuno from University of Tsukuba.