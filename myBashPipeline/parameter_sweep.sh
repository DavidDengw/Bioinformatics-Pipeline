#!/bin/bash

cd ./tiny-test-data/wgs ||
#Parameter Sweep on -k Minimum seed length, comparing average mapping quality.
for i in 10 20 30; do bwa mem -k $i ../genomes/Hsapiens/hg19/bwa/hg19.fa mt_1.fq.gz mt_2.fq.gz > mt_k$i.sam; done

#Filter out header lines that start with '@', only keep parameter values and mapping quality.
for f in mt_k*; do echo $f; grep -v '^@' $f | awk '{ sum += $5; n++ } END { if (n > 0) print sum / n; }'; done | sed 's/mt_k//' | sed 's/.sam//'


