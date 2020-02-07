#!/bin/sh
if test "$1" = ""
then echo "Sample list?"
     read list
else list=$1
fi

for sample in `cat $list`
do hisat2 -x ~/Genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2Index/genome -1 $sample'_R1.fastq' -2 $sample'_R2.fastq' -S Mapping_default_parameters_$sample'.sam'
done
