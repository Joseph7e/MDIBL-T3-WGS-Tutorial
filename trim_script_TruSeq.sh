#!/bin/bash

FORWARD=$1
REVERSE=$2

trimmomatic PE -threads 72 $FORWARD $REVERSE\
    paired_forward.fastq.gz unpaired_forward.fastq.gz\
    paired_reverse.fastq.gz unpaired_reverse.fastq.gz\
    ILLUMINACLIP:/usr/share/trimmomatic/TruSeq2-PE.fa:2:30:10\
    SLIDINGWINDOW:4:15 MINLEN:36
