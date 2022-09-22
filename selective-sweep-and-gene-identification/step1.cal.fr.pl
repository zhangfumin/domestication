#!/usr/bin/perl -w

#instruction: Three parameters must be provided on the command line, sample, vcf, fr.
#sample: a file contains the informaton of each individual sample and its population. Every line includes names of an individual sample and its population, whiche were seperated by a tab character. If samples are used in next analysis, the sign '+' were add in front of teir name. For example, +ind_CHN_MC18	cultivated 
#vcf: an vcf file contains SNPs of all samples
#fr: an Output fule including SNP allele frequency in ecah population

use strict;

my $sample = shift;
my $vcf    = shift;   
my $fr     = shift;

my ($i,$j);
my @tm=();

my @pop_name=("cultivated","wild");
my %pop_pt=();
foreach(@pop_name){$pop_pt{$_}=1;}
my %pop_sam=();

open (R, $sample) || die "reading $sample error..\n";
while(<R>)
  { if($_ !~ /^\+/)  { next;}
    chomp $_;
    $_ =~ s/^\+//;
    @tm=split(/\t/,$_);
    if(exists($pop_pt{$tm[1]})){push(@{$pop_sam{$tm[1]}}, $tm[0]); push(@{$pop_sam{"all"}}, $tm[0]);}
  }
close(R);  

my %sam_pos=();
open(R,$vcf) || die "Reading $vcf failed ..\n";
open(W,">".$fr) || die "Writing $vcf failed ..\n";
while(<R>)
     { chomp $_;
       if($_ =~ /^\#\#/){next;}
       if($_ =~ /^\#CHROM/) 
         { @tm = split(/\t/, $_);
           printf W "%s\t%s\t%s\t%s\t%s", $tm[0],$tm[1],$tm[3],$tm[4],"ALL";
           foreach(@pop_name) {printf W "\t%s", $_;}
           print W "\n";
           for $i(9..scalar(@tm)-1) {$sam_pos{$tm[$i]}=$i;}
           next;
          }
       @tm = split(/\t/, $_);
       $tm[0]=~ s/build05r1.fasta//;
       if($tm[4] =~ /\*/) {next;}
       my @allele=split(/,/,$tm[4]);
       my $aln=scalar(@allele);
       printf W "%s\t%s\t%s\t%s", $tm[0], $tm[1], $tm[3], $tm[4];
       my @seg=();
       foreach(@{$pop_sam{'all'}}) { push(@seg, $tm[$sam_pos{$_}]);} 
       my @fr=cal_fr($aln, \@seg);
       my $frn=$fr[0];
       for $i(1..scalar(@fr)-1){$frn .= ':'.$fr[$i];}
       printf W "\t%s", $frn;
       foreach(@pop_name)
              { @seg=();
                foreach(@{$pop_sam{$_}}) {push(@seg, $tm[$sam_pos{$_}]);}
                @fr=cal_fr($aln, \@seg);
                $frn=$fr[0];
                for $i(1..scalar(@fr)-1){$frn .= ':'.$fr[$i];}
                printf W "\t%s", $frn;
              }  
       print W "\n"; 
       }                                              
close(R);
close(W);
         
sub   cal_fr  
       { my ($an, $seg)=@_;
         my %ale=();
         my $i=0;
         for $i(0..$an) {$ale{$i}=0;}
         foreach (@$seg)
                 { if( $_ =~ /^\./) {next;} 
                   my @tk=split(/:/,$_);
                   my @tn=();
                   if($tk[0] =~ /\//){@tn=split(/\//, $tk[0]);}
                   else              {@tn=split(/\|/, $tk[0]);}                           
                   $ale{$tn[0]}++; 
                   $ale{$tn[1]}++;
                 }
         my @sfr=();
         for  $i(0..$an) { push (@sfr, $ale{$i}); }
         return @sfr;
        }
                 
                         
             


