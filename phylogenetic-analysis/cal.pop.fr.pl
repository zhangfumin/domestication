#!/usr/bin/perl -w

#instruction: Three parameters must be provided on the command line, saminfo, vcf_file, fre_file.
#sampinfo: a file contains the informaton of each individual sample and its population. Every line includes names of an individual sample and its population, whiche were seperated by a tab character. If one individual sample is not used in next analysis, xxx can be used as its population name.  
#vcf: an vcf file contains SNPs of all samples
#fre_file: an Output fule including SNP allele frequency in ecah population                     

use strict;

my $sam_info = shift;
my $vcf_file = shift;   
my $fre_file = shift;

my ($i,$j);
my @tm=();

my %pop_idx = ();
my @pop_nam = ();
my %sam_idx = ();

open(R,$sam_info) || die "reading $sam_info failed ..\n";
while(<R>)
  { chomp $_;
    @tm=split(/\t/,$_);
    if( $tm[1] ne 'xxx' && length($tm[1])>0 )
      {$sam_idx{$tm[0]}=$tm[1];
       if(!exists($pop_idx{$tm[1]}))
         {$pop_idx{$tm[1]}=1; push(@pop_nam,$tm[1]);}
      }
  }    
close(R);

my %pop_sam=();
   @{$pop_sam{'ALL'}}=();
foreach(@pop_nam){@{$pop_sam{$_}}=();}   
my %sam_pos=();
open(R,$vcf_file) || die "Reading $vcf_file failed ..\n";
open(W,">".$fre_file) || die "Writing $fre_file failed ..\n";
while(<R>)
     { chomp $_;
       if($_ =~ /^\#\#/){next;}
       if($_ =~ /^\#CHROM/) 
         { @tm = split(/\t/, $_);
           printf W "%s\t%s\t%s\t%s\t%s", $tm[0],$tm[1],$tm[3],$tm[4],"ALL";
           foreach(@pop_nam) {printf W "\t%s", $_;}
           print W "\n";
           for $i(9..scalar(@tm)-1) 
                 { $sam_pos{$tm[$i]}=$i;
                   if(exists($sam_idx{$tm[$i]}))
                     { push(@{$pop_sam{'ALL'}},$tm[$i]);
                       push(@{$pop_sam{$sam_idx{$tm[$i]}}},$tm[$i]);
                     }
                 }    
           next;
         }
       @tm = split(/\t/, $_);
       if($tm[4] =~ /\*/) {next;}
       my @allele=split(/,/,$tm[4]);
       my $aln=scalar(@allele);
       printf W "%s\t%s\t%s\t%s", $tm[0], $tm[1], $tm[3], $tm[4];
       my @seg=();
       foreach(@{$pop_sam{'ALL'}}) { push(@seg, $tm[$sam_pos{$_}]);} 
       my @fr=cal_fr($aln, \@seg);
       my $frn=$fr[0];
       for $i(1..scalar(@fr)-1){$frn .= ':'.$fr[$i];}
       printf W "\t%s", $frn;
       foreach(@pop_nam)
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
                 
                         
             


