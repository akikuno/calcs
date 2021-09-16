#!/bin/sh

ref=tests/random_100bp.fa
que=tests/deletion/del.fq

minimap2 --cs "$ref" "$que" >tests/deletion/del_cs.paf
minimap2 --cs=long "$ref" "$que" >tests/deletion/del_cslong.paf

minimap2 -ax map-ont "$ref" "$que" >tests/deletion/del.sam
minimap2 -ax map-ont --cs "$ref" "$que" >tests/deletion/del_cs.sam
minimap2 -ax map-ont --cs=long "$ref" --cs "$que" >tests/deletion/del_cslong.sam
