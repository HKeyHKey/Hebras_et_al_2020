#!/usr/bin/perl

if ($ARGV[1] eq '')
{
    print "Please enter sample name, then name of BED file containing selected read pairs (e.g., ./Module_converts_quality_scores.pl OK6_GTGCCA_L008 Amplicon_reads_from_Mapping_default_parameters_OK6_GTGCCA_L008.bed).\n";
}
else
{
    $s=0;
    foreach $char ('!','"','#','$','%','&','\'','(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J')
    {
        $score{$char}=$s;
        ++$s;
    }

    open(BED,$ARGV[1]);
    while (<BED>)
    {
        chomp;
        @array=split('\t',$_);
        if ($array[3]=~/\/2$/)
        {
            $ID=$array[3];
            $ID=~s/\/2$//;
            push(@select,$ID);
        }
    }
    close(BED);

    open(IN,"Mapping_default_parameters_".$ARGV[0].".sam");
    open(OUT,">Quality_in_".$ARGV[0].".dat");
    print OUT "Read_ID";
    for ($i=1;$i<=150;++$i)
    {
        print OUT " nt_$i";
    }
    for ($i=1;$i<=150;++$i)
    {
        print OUT " q_$i";
    }
    print OUT "\n";


    while (<IN>)
    {
        chomp;
        if (/^\@PG\t/)
        {
            $read=1;
        }
        else
        {
            if ($read)
            {
                @array=split('\t',$_);
                ($ID,$seq,$qual)=($array[0],$array[9],$array[10]);
                @seq_array=split('',$seq);
                @qual_array=split('',$qual);
                if ($ID eq $select[$select_index]) # Here I am taking advantage of the fact that reads are sorted in the same order in the BED and SAM files (this wouldn't work otherwise!)
                {
                    if ($ID eq $memo) #So we only print sequence and quality for "R2" reads (not their "R1" mate in the read pair)
                    {
                        ++$select_index;
                        print OUT "$ID";
                        foreach $nt (@seq_array)
                        {
                            print OUT " $nt";
                        }
                        foreach $q (@qual_array)
                        {
                            print OUT " $score{$q}";
                        }
                        print OUT "\n";
                    }
                    else
                    {
                        $memo=$ID;
                    }
                }
            }
        }
    }
    close(IN);
    close(OUT);
}
