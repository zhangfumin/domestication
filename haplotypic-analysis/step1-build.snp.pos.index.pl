#!/usr/bin/perl -w

use strict;

my $vcf = shift; # the VCF file contains the data of whole-genome SNPs

my %chr_pos = ();
my @chr_key=();
my @tm=();
my $pos_marker= 0; 
my $pos_value = 0;
open(R,$vcf) || die "reading $vcf failed .. \n";

while(<R>)
     { if( $_ =~ /\#/) {$pos_value = tell(R); next;}
       @tm = split(/\t/,$_);
       $pos_marker = int($tm[1]/1000000);
       if( !exists($chr_pos{$tm[0]."-".$pos_marker}) )
         { $chr_pos{$tm[0]."-".$pos_marker} = $pos_value;
           push(@chr_key, $tm[0]."-".$pos_marker);
          }
       $pos_value = tell(R);  
      }
open(W,">chrom.snp.pos");
foreach (@chr_key) {printf W "%s\t%s\n", $_, $chr_pos{$_};}
close(W);   
