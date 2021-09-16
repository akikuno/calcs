#!/bin/sh

ref=tests/random_100bp.fa
que=tests/subindel/subindel.fq

rm tests/subindel/*.paf tests/subindel/*.sam

# PAF
minimap2 --cs "$ref" "$que" >tests/subindel/subindel_cs.paf
minimap2 --cs=long "$ref" "$que" >tests/subindel/subindel_cslong.paf

# SAM
minimap2 -ax map-ont "$ref" "$que" >tests/subindel/subindel.sam
minimap2 -ax map-ont --cs "$ref" "$que" >tests/subindel/subindel_cs.sam
minimap2 -ax map-ont --cs=long "$ref" "$que" >tests/subindel/subindel_cslong.sam

# Check CS tag
cat tests/subindel/subindel_cs.sam | awk '$1 !~ "@" {print $(NF-1)}'
cat tests/subindel/subindel_cslong.sam | awk '$1 !~ "@" {print $(NF-1)}'
