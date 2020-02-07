#!/bin/sh

if test "$1" = ""
then echo "Sample list?"
     read list
else list=$1
fi

for sample in `cat $list`
do Rscript R_editing_combinations $sample
done
