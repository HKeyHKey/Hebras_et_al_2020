#!/bin/sh

for f in `ls Quality_in_$1''*.dat`
do name=`echo $f | sed -e 's|Quality_in_||' -e 's|\.dat||'`
   Rscript R_measures_editing $name
done
