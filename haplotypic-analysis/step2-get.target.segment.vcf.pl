#!/usr/bin/perl -w

use strict;

my  $gene     =  shift;   #  a list containing inforamtion of taget genes, inculding the fields of name, chromosome, start position, end position. pleas note that the fields were seperated by delimiter "TAB"  
my  $up_down  =  shift;   #  how many MB of up- and down-stream of the target genes were used for the next analysis   
my  $segdir   =  shift;   #  a folder for saving results. if the folder doesn't exist, it will be established.
my  $vcf_snp  =  shift;   #  the VCF file contans the data of  whole-genome SNPs

my $time0 = time();

if(!opendir(SEG, $segdir)){ mkdir($segdir,oct(750)) || die "Establishing the folder $segdir failed .. /n"; }
else                      { closedir(SEG); }

my   @tm =();
my   ($i, $j);
my   @gen_mat = ();
open (R,$gene) || die "Couldn't open $gene for reading .. \n";
$i=0;
while(<R>)
     { chomp $_;
       $_ =~  s/ //g;
       @tm=split(/\t/,$_);
       push(@{$gen_mat[$i]},@tm);
       $i++;
     }
close(R);

my @sam=();
open(R,"201110.sample.grp") || die "Couldn't open 200928.1493.samples.grp for reading ..\n";
while ( <R> )
      { chomp $_;
        @tm=split(/\t/, $_);
        push(@sam, $tm[0]);
       }
close(R);

my %chr_snp_pos = ();
open(R, "chrom.snp.pos") || die "Couldn't open chrom.pos for reading .. ";
while ( <R> )
      { chomp $_;
        @tm=split(/\t/,$_);
        $chr_snp_pos{$tm[0]} = $tm[1];
      }
close(R);

open(R, $vcf_snp)  ||  die "Couldn't open  $vcf_snp for reading ..\n ";

my @headseg=();
my @headline=();
my %sam_snp_pos = ();
while( <R> )   
    { if ($_ =~ /^\#\#/ )       { push(@headseg, $_); }
      if ($_ =~ /^\#CHROM/ )    { chomp $_; @tm=split(/\t/,  $_);   last; } 
    } 
for $i(0..8) { push(@headline, shift(@tm)); }
$i=9;  foreach(@tm) {$sam_snp_pos{$_}=$i; $i++;}

foreach my $gm( @gen_mat) { open(W,">".$segdir."/". $gm->[0]. ".segment.vcf");
                            foreach(@headseg) {print W $_;}
                            printf W "%s", $headline[0];
                            for $i(1..8)  { printf W "\t%s",$headline[$i]; }
                            foreach(@sam) { printf W "\t%s", $_; }
                            print W "\n";
                            my $ss =$gm->[2]-$up_down*1000000;
                            my $es =$gm->[3]+$up_down*1000000;;
                            if($ss<0) {$ss=0;}
                            my $pos_marker=int($ss/1000000);
                            seek(R, $chr_snp_pos{ $gm->[1].'-'.$pos_marker},0);
                            my @pseg=();
                            while (<R>)
                                  {  chomp $_;
                                     @tm=split(/\t/, $_);
                                     if ($tm[0] eq $gm->[1]  &&  $ss <= $tm[1] && $tm[1] <= $es) {push(@pseg, $_);}
                                     if ($tm[1] > $es || $tm[0] ne $gm->[1])  { last; }
                                  } 
                            foreach (@pseg)
                                 {  @tm = split (/\t/, $_);
                                    printf W "%s", $tm[0];
                                    for $i(1..8) { printf W "\t%s",$tm[$i]; }
                                    foreach (@sam)  { printf W "\t%s",$tm[$sam_snp_pos{$_}]; }
                                    print W "\n";
                                 } 
                            close(W); 
                        }
close(R);
my $time1 = time();
my $tsec=$time1-$time0;
my $tmin=int($tsec/60);
my $sec=$tsec%60;
printf "Elapsed time : %s min %s sec\n", $tmin, $sec; 

                                          
