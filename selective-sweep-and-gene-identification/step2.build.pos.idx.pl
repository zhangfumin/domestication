#!/usr/bin/perl -w

use strict;

my $fr   = shift;
my $idx  = shift;

my $ccr="chr0";
my %chr_pos = ();
my @chr_key=();
my @tm=();
my $pos_marker= 0; 
my $pos_value = 0;
open(R,$fr);
while(<R>)
     { if( $_ =~ /\#/) {$pos_value = tell(R); next;}
       @tm = split(/\t/,$_);
       if($tm[0] ne $ccr ) 
         {$ccr=$tm[0]; 
          $chr_pos{$tm[0]."-0"}=$pos_value; 
          push(@chr_key, $tm[0]."-0");
         }
       $pos_marker = int($tm[1]/10000);
       if(!exists($chr_pos{$tm[0]."-".$pos_marker})) 
         { $chr_pos{$tm[0]."-".$pos_marker} = $pos_value;
           push(@chr_key, $tm[0]."-".$pos_marker);   
         }
       $pos_value = tell(R);  
      }
open(W,">".$idx);
foreach (@chr_key) {printf W "%s\t%s\n", $_, $chr_pos{$_};}
close(W);   
