#!/bin/sh

#Conversion of the SAM files into BAM (simultaneously selecting reads mapped in proper pairs using option "-f 2" and excluding non-primary alignments using "-F 256") then BED, then intersection with amplicon BED file:
for file in `ls Mapping*.sam`;do samtools view -f 2 -F 256 -Sb $file > `echo $file | sed 's|\.sam|.bam|'`;done
for file in `ls Mapping*.bam`;do bedtools bamtobed -i $file > `echo $file | sed 's|\.bam|.bed|'`;done
for file in `ls Mapping*.bed`;do bedtools intersect -wa -a $file -b Amplicon.bed > Amplicon_reads_from_$file;done

