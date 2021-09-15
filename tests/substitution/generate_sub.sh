#!/bin/sh

ref=tests/random_100bp.fa
que=tests/substitution/sub.fq

minimap2 --cs "$ref" "$que" >tests/substitution/sub_cs.paf
minimap2 --cs=long "$ref" "$que" >tests/substitution/sub_cslong.paf

minimap2 -ax map-ont --cs "$ref" "$que" >tests/substitution/sub_cs.sam
minimap2 -ax map-ont --cs=long "$ref" --cs "$que" >tests/substitution/sub_cslong.paf
