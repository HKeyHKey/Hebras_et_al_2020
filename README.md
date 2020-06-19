# Hebras_et_al_2020
Scripts and intermediary data for analyses published in Figures 3, 4, 5, 6 and Supplementary Figure 2 of Hebras et al. (2020).

Two series of samples were prepared and sequenced 9 months apart (see file 'Dissected_area_nomenclature_in_Figure3.txt' for the assignation of the samples shown in Figure 3 into these two series).

## 1. Fig. 3: Analysis of the first series of samples (SRA accession #PRJNA603264):

``mkdir First_series;cd First_series``

(then copy every file in this directory, except files annotated '_series_2')

``cp Dissected_areas_series_1.txt Dissected_areas.txt``

### 1.1. Aligning on the mm10 genome:

``ls *.fastq | sed 's|_R[12]\.fastq||' | sort | uniq > sample_list
for i in `seq 1 8`;do j=`echo "("$i"-1)*26+1" | bc`;tail -n +$j sample_list | head -26 > sample_list_$i;done
for i in `seq 1 8`;do ./Script_align.sh sample_list_$i > nohup_$i'.out';sleep 3;done``

### 1.2. Conversion of the SAM files into BAM (simultaneously selecting reads mapped in proper pairs using option "-f 2" and excluding non-primary alignments using "-F 256") then BED, then intersection with amplicon BED file:

``./Script_intersect.sh``

### 1.3. Collecting library statistics:

``wc -l *.fastq > File_size_series_1.dat;echo "Sample Number_of_read_pairs Number_of_nonconcordant_pairs Number_of_uniquely_matching_concordant_pairs Number_of_multiply_matching_concordant_pairs Number_of_amplicon_matching_concordant_pairs" > Library_statistics_series_1.dat;for sample in `ls *.sam | sed -e 's|Mapping_default_parameters_||' -e 's|\.sam$||'`;do chunk=`grep -nw $sample sample_list_* | sed -e 's|^sample_list_||' -e 's|:.*||'`;n=`grep -nw $sample sample_list_* | sed -e 's|.*:\(.*\):.*|\1|'`;grep -A 4 -m $n 'reads; of these:' nohup_$chunk'.out' | tail -5 > tmp;total=`head -1 tmp | awk '{print $1}'`;nonconcordant=`head -3 tmp | tail -1 | awk '{print $1}'`;concordant1=`head -4 tmp | tail -1 | awk '{print $1}'`;concordantmult=`tail -1 tmp | awk '{print $1}'`;amplicon2=`cat Amplicon_reads_from_Mapping_default_parameters_$sample'.bed' | wc -l`;amplicon=`echo $amplicon2"/2" | bc`;echo $sample $total $nonconcordant $concordant1 $concordantmult $amplicon;done >> Library_statistics_series_1.dat``

### 1.4. Selection of high-quality nucleotides (quality score >= 20) and evaluation of editing frequencies on these reads:

``ls Amplicon_reads_from_Mapping_default_parameters_* > amplicon_list;for i in `seq 1 8`;do j=`echo "("$i"-1)*26+1" | bc`;tail -n +$j amplicon_list | head -26 > amplicon_list_$i;done;for i in `seq 1 8`;do ./Script_quality.sh amplicon_list_$i > nohup_amplicon_$i'.out';sleep 3;done``

(can be launched in parallel on 8 processors instead of launching them sequentially with the "for" loop)

### 1.5. Extracting frequencies of A-I edition:

``for f in `ls Quality_in_i*.dat`;do name=`echo $f | sed 's|Quality_in_i|Quality_in_I|'`;mv $f $name;done #To homogeneize case``

``for code in O H I S C T X R W;do ./Script_measure.sh $code > nohup_$code'.txt';sleep 3;done``

(can be launched in parallel on 9 processors instead of launching them sequentially with the "for" loop)

``for f in `ls Editing_data_*.csv`;do sed -e 's|^"[0-9]*",||' -e "1 s|.*|Position in 'genomic_sequence.fa',Number of reads,Number of reads with A,Number of reads with C, Number of reads with G, Number of reads with T|" $f > tmp_$f;mv -f tmp_$f $f;done``

Resulting files: 'Editing_data_*.csv', stored in archive 'Series_1_editing_data.tar.bz2'; these files contain the number of reads with A, C, G or T at every position (between nt 63 and 130 of 'genomic_sequence.fa') where the genome contains an A.

``for area in `tail -n +2 Dissected_areas_series_1.txt | awk '{print $1}'`;do Rscript R_editing_stats $area;done;for f in `ls Values_for_plotting_editing_frequency_*.csv`;do sed -e 's|m_wt|mean (wt)|' -e 's|se_wt|st. error (wt)|' -e 's|m_mut|mean (mutant)|' -e 's|se_mut|st. error (mutant)|' $f > tmp;mv -f tmp $f;done``

Resulting files: 'Values_for_plotting_editing_frequency_*.csv', stored in archive 'Series_1_editing_frequency_site_by_site.tar.bz2'. Values for olfactory bulb (sample code: O), hypothalamus (sample code: H), hippocampus (sample code: I), striatum (sample code: S), cerebellum (sample code: C), brainstem (sample code: T), prefrontal cortex (sample code: X) and raphe nuclei (sample code: R) were used to plot Figure 3B.


## 2. Fig. 3 and S2: Analysis of the second series of samples (SRA accession #PRJNA603261):

``mkdir ../Second_series;cd ../Second_series``

(then copy every file in this directory, except files annotated '_series_1')

``cp Dissected_areas_series_2.txt Dissected_areas.txt``


### 2.1. Aligning on the mm10 genome:

``ls *.fastq | sed 's|_R[12]\.fastq||' | sort | uniq > sample_list;for i in `seq 1 8`;do j=`echo "("$i"-1)*24+1" | bc`;tail -n +$j sample_list | head -24 > sample_list_$i;done;for i in `seq 1 8`;do ./Script_align.sh sample_list_$i > nohup_$i'.out';sleep 3;done``

### 2.2. Conversion of the SAM files into BAM (simultaneously selecting reads mapped in proper pairs using option "-f 2" and excluding non-primary alignments using "-F 256") then BED, then intersection with amplicon BED file:

``./Script_intersect.sh``

### 2.3. Collecting library statistics:

``wc -l *.fastq > File_size_series_2.dat;echo "Sample Number_of_read_pairs Number_of_nonconcordant_pairs Number_of_uniquely_matching_concordant_pairs Number_of_multiply_matching_concordant_pairs Number_of_amplicon_matching_concordant_pairs" > Library_statistics_series_2.dat;for sample in `ls *.sam | sed -e 's|Mapping_default_parameters_||' -e 's|\.sam$||'`;do chunk=`grep -nw $sample sample_list_* | sed -e 's|^sample_list_||' -e 's|:.*||'`;n=`grep -nw $sample sample_list_* | sed -e 's|.*:\(.*\):.*|\1|'`;grep -A 4 -m $n 'reads; of these:' nohup_$chunk'.out' | tail -5 > tmp;total=`head -1 tmp | awk '{print $1}'`;nonconcordant=`head -3 tmp | tail -1 | awk '{print $1}'`;concordant1=`head -4 tmp | tail -1 | awk '{print $1}'`;concordantmult=`tail -1 tmp | awk '{print $1}'`;amplicon2=`cat Amplicon_reads_from_Mapping_default_parameters_$sample'.bed' | wc -l`;amplicon=`echo $amplicon2"/2" | bc`;echo $sample $total $nonconcordant $concordant1 $concordantmult $amplicon;done >> Library_statistics_series_2.dat``

### 2.4. Selection of high-quality nucleotides (quality score >= 20) and evaluation of editing frequencies on these reads:

``ls Amplicon_reads_from_Mapping_default_parameters_* > amplicon_list;for i in `seq 1 8`;do j=`echo "("$i"-1)*24+1" | bc`;tail -n +$j amplicon_list | head -24 > amplicon_list_$i;done;for i in `seq 1 8`;do ./Script_quality.sh amplicon_list_$i > nohup_amplicon_$i'.out';sleep 3;done``

(can be launched in parallel on 8 processors instead of launching them sequentially with the "for" loop)

### 2.5. Extracting frequencies of A-I edition:

``for code in d e h m M N P w;do ./Script_measure.sh $code > nohup_$code'.txt';sleep 3;done``

(can be launched in parallel on 8 processors instead of launching them sequentially with the "for" loop)

``for f in `ls Editing_data_e16k*`;do name=`echo $f | sed 's|data_e16k\([0-9]*\)_|data_e16\1k_|'`;mv $f $name;done #To homogeneize nomenclature``

``for f in `ls Editing_data_*.csv`;do sed -e 's|^"[0-9]*",||' -e "1 s|.*|Position in 'genomic_sequence.fa',Number of reads,Number of reads with A,Number of reads with C, Number of reads with G, Number of reads with T|" $f > tmp_$f;mv -f tmp_$f $f;done``

Resulting files: 'Editing_data_*.csv', stored in archive 'Series_2_editing_data.tar.bz2'; these files contain the number of reads with A, C, G or T at every position (between nt 63 and 130 of 'genomic_sequence.fa') where the genome contains an A

``for area in `tail -n +2 Dissected_areas_series_2.txt | awk '{print $1}'`;do Rscript R_editing_stats $area;done;for f in `ls Values_for_plotting_editing_frequency_*.csv`;do sed -e 's|m_wt|mean (wt)|' -e 's|se_wt|st. error (wt)|' -e 's|m_mut|mean (mutant)|' -e 's|se_mut|st. error (mutant)|' $f > tmp;mv -f tmp $f;done``

Resulting files: 'Values_for_plotting_editing_frequency_*.csv', stored in archive 'Series_2_editing_frequency_site_by_site.tar.bz2'. Values for whole brain (sample code: wb) and spinal chord (sample code: ME) were used to plot Figure 3B.


## 3. Analysis of combinations of site edition (Figure 3C):

### 3.1. File preparation for the first series of samples (SRA accession #PRJNA603264):

``cd ../First_series;for i in `seq 1 8`;do ./Script_editing_combination_measure.sh sample_list_$i > nohup_ter_$i'.out';sleep 3;done``

(can be launched in parallel on 8 processors instead of launching them sequentially with the "for" loop)

``for area in `tail -n +2 Dissected_areas_series_1.txt | awk '{print $1}'`;do for f in `ls Editing_combinations_in_$area''*`;do name=`echo $f | sed -e 's|Editing_combinations_in_||' -e 's|_.*||'`;echo "Combination Number_in_"$name > tmp_$f;sed '/"/ {
N
s|^\[1\] *"\(.*\)"\n\[1\] *\(.*\)|\1 \2|
}' $f >> tmp_$f;done``

### 3.2. File preparation for the second series of samples (SRA accession #PRJNA603261):

``cd ../Second_series;for i in `seq 1 8`;do ./Script_editing_combination_measure.sh sample_list_$i > nohup_ter_$i'.out';sleep 3;done``

(can be launched in parallel on 8 processors instead of launching them sequentially with the "for" loop)

``for area in `tail -n +2 Dissected_areas_series_1.txt | awk '{print $1}'`;do for f in `ls Editing_combinations_in_$area''*`;do name=`echo $f | sed -e 's|Editing_combinations_in_||' -e 's|_.*||'`;echo "Combination Number_in_"$name > tmp_$f;sed '/"/ {
N
s|^\[1\] *"\(.*\)"\n\[1\] *\(.*\)|\1 \2|
}' $f >> tmp_$f;done``

### 3.3 Joint analysis of the 10 samples shown on Figure 3 (adult whole brain and 9 adult brain areas, taken from both series):

Nomenclature homogenization: first series:

``cd ..;for code in O H I S C T X R;do for f in `ls First_series/tmp_Editing_combinations_in_$code''[0-9kK]*| grep -v '_series/tmp_Editing_combinations_in_'$code'[0-9]*[kK]'`;do replicate=`echo $f | sed -e 's|First_series/tmp_Editing_combinations_in_'$code'||' -e 's|_.*||'`;sed '1 s|^|#|' $f > Combinations_$code'_wt_'$replicate'.dat';done;for f in `ls First_series/tmp_Editing_combinations_in_$code''[0-9kK]*| grep '_series/tmp_Editing_combinations_in_'$code'[0-9]*[kK]'`;do replicate=`echo $f | sed -e 's|First_series/tmp_Editing_combinations_in_'$code'||' -e 's|_.*||' -e 's|[Kk]||'`;sed '1 s|^|#|' $f > Combinations_$code'_mut_'$replicate'.dat';done;done``

second series:

``for code in ME wb e16 e16h P1H e16c e16i hypo del32H del32P1H mCPP Nacl;do for f in `ls Second_series/tmp_Editing_combinations_in_$code''[0-9kK]* | grep -v '_series/tmp_Editing_combinations_in_'$code'[0-9]*[kK]'`;do replicate=`echo $f | sed -e 's|Second_series/tmp_Editing_combinations_in_'$code'||' -e 's|_.*||'`;sed '1 s|^|#|' $f > Combinations_$code'_wt_'$replicate'.dat';done;for f in `ls Second_series/tmp_Editing_combinations_in_$code''[0-9kK]* | grep '_series/tmp_Editing_combinations_in_'$code'[0-9]*[kK]'`;do replicate=`echo $f | sed -e 's|Second_series/tmp_Editing_combinations_in_'$code'||' -e 's|_.*||' -e 's|[Kk]||'`;sed '1 s|^|#|' $f > Combinations_$code'_mut_'$replicate'.dat';done;done``

Files 'Combinations_*.dat' are stored in archive 'Editing_combination_counts_replicate_per_replicate.tar.bz2'.

``R CMD BATCH R_commands_combination_analysis``

Resulting file: 'Combination_comparisons.csv', used in Figure 3C.

## 4. ANOVA for the effect of development, brain area, editing site identity, and genotype (for Fig. S2):

### 4.1. Nomenclature homogenization:

``for code in O H I S C T X R ME wb e16 e16h P1H e16c e16i hypo del32H del32P1H mCPP Nacl;do for f in `ls Editing_data_$code''[0-9kK]* | grep -v Editing_data_$code''[0-9]*[kK]`;do replicate=`echo $f | sed -e 's|^Editing_data_'$code'||' -e 's|_.*||' -e 's|[kK]||'`;mv $f Site_by_site_$code'_wt_'$replicate'.csv';done;done;for code in O H I S C T X R ME wb e16 e16h P1H e16c e16i hypo del32H del32P1H mCPP Nacl;do for f in `ls Editing_data_$code''[0-9kK]* | grep Editing_data_$code''[0-9]*[kK]`;do replicate=`echo $f | sed -e 's|^Editing_data_'$code'||' -e 's|_.*||' -e 's|[kK]||'`;mv $f Site_by_site_$code'_mut_'$replicate'.csv';done;done``

Resulting files ('Site_by_site_*.csv') can be found in archive 'Data_for_Fig_S2.tar.bz2'.

### 4.2 ANOVA:

``R CMD BATCH R_editing_in_development``

## 5. ANOVA for Fig. 4 (panels A to F), Fig. 5 (panels C to L) and Fig. 6 (panels A to G):

Categorical variables (multivariate histograms in Figures 4, 5 and 6):

``for f in `ls Fig*.csv | grep -v Fig[56]C | grep -v Fig6A`;do Rscript R_commands_ANOVA_histograms $f;done``

Time course analyses with non-linear time-dependency (Fig. 5C and 6C):

``for f in Fig5C.csv Fig6C.csv;do Rscript R_commands_ANOVA_time_series_non_linear $f;done``

Time course analysis with linear time-dependency (Fig. 6A):

``Rscript R_commands_ANOVA_time_series_linear Fig6A.csv``
