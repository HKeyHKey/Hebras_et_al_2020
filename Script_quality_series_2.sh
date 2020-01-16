#!/bin/sh

for f in `cat $1`;do name=`echo $f | sed -e 's|Amplicon_reads_from_Mapping_default_parameters_||' -e 's|\.bed||'`;./Module_converts_quality_scores.pl $name $f;Rscript R_display_quality $name;done

