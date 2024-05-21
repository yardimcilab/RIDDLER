#!/bin/sh -ex
#Calculate %GC for fixed size bins in assembly fasta
#Output bed contains other features, GC% is 5th column

resolution=$1
assembly=$2
binBed=$3

bedtools nuc -fi $assembly -bed $binBed > GCcontent.$resolution.nuc.bed

