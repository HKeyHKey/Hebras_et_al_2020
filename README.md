# Hebras_et_al_2020
Scripts and intermediary data for analyses published in Figure 3 of Hebras et al. (2020).

Two series of samples were prepared and sequenced 9 months apart (see file 'Dissected_area_nomenclature.txt' for the assignation of the samples shown in Figure 3 into these two series).

## 1. Analysis of the first series of samples (SRA accession #???):

### 1.1. Aligning on the mm10 genome:

``ls *.fastq | sed 's|_R[12]\.fastq||' | sort | uniq > sample_list
for i in `seq 1 8`;do j=`echo "("$i"-1)*26+1" | bc`;tail -n +$j sample_list | head -26 > sample_list_$i;done
for i in `seq 1 8`;do ./Script_align_series_1.sh sample_list_$i > nohup_$i'.out';sleep 3;done``

### 1.2. Conversion of the SAM files into BAM (simultaneously selecting reads mapped in proper pairs using option "-f 2" and excluding non-primary alignments using "-F 256") then BED, then intersection with amplicon BED file:

``./Script_intersect_series_1.sh``

### 1.3. Collecting library statistics:

``wc -l *.fastq > File_size.dat;echo "Sample Number_of_read_pairs Number_of_nonconcordant_pairs Number_of_uniquely_matching_concordant_pairs Number_of_multiply_matching_concordant_pairs Number_of_amplicon_matching_concordant_pairs" > Library_statistics.dat;for sample in `ls *.sam | sed -e 's|Mapping_default_parameters_||' -e 's|\.sam$||'`;do chunk=`grep -nw $sample sample_list_* | sed -e 's|^sample_list_||' -e 's|:.*||'`;n=`grep -nw $sample sample_list_* | sed -e 's|.*:\(.*\):.*|\1|'`;grep -A 4 -m $n 'reads; of these:' nohup_$chunk'.out' | tail -5 > tmp;total=`head -1 tmp | awk '{print $1}'`;nonconcordant=`head -3 tmp | tail -1 | awk '{print $1}'`;concordant1=`head -4 tmp | tail -1 | awk '{print $1}'`;concordantmult=`tail -1 tmp | awk '{print $1}'`;amplicon2=`cat Amplicon_reads_from_Mapping_default_parameters_$sample'.bed' | wc -l`;amplicon=`echo $amplicon2"/2" | bc`;echo $sample $total $nonconcordant $concordant1 $concordantmult $amplicon;done >> Library_statistics.dat``

### 1.4. Selection of high-quality nucleotides (quality score >= 20) and evaluation of editing frequencies on these reads:

``./Script_quality_series_1.sh # launched on sub-lists of all the 'Amplicon_reads_from_Mapping_default_parameters_*.bed' files; to extract quality scores and convert from ASCII into numeric, and to plot quality scores along "R2" reads in amplicon-matching concordant pairs (these graphs show that the chosen cutoff on quality score, 20, selects most reads at each of the 5 editing sites)``

### 1.5. Extracting frequencies of A-I edition:

``for f in `ls Quality_in_i*.dat`;do name=`echo $f | sed 's|Quality_in_i|Quality_in_I|'`;mv $f $name;done #To homogeneize case;for code in O H I S C T X R W;do ./Script_measure_series_1.sh O  > nohup_O.txt;sleep 3;done``

Resulting files: 'Editing_data_*.csv'. stored in archive 'Series_1_editing_data.tar.bz2'. Values for olfactory bulb (sample code: O), hypothalamus (sample code: H), hippocampus (sample code: I), striatum (sample code: S), cerebellum (sample code: C), brainstem (sample code: T), prefrontal cortex (sample code: X) and raphe nuclei (sample code: R) were used to plot Figure 3B.
