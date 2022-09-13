#!/usr/bin/perl -w

#In this script, we used the progam neighbpor of the sofware phylip and made a little modification in source C code to remove intactive interface

use strict;

my $in        = shift; #file including pi, dxy and fst of gene or window calculated by cal_pop_parameter.pl
my $dxy_tree  = shift; #outfile 
my $fst_tree  = shift; #outfile

my @sample=qw(R1 JAP IND N2 N1 R2 R2B); #group name
my $sn=scalar(@sample);
my $dxy_sta=5+$sn;
my $dxy_end=$dxy_sta+$sn*($sn-1)/2-1;
my $fst_sta=$dxy_end+1;
my $fst_end=$fst_sta+$sn*($sn-1)/2-1;

open(WDT,">".$dxy_tree);
open(WFT,">".$fst_tree);

my($i,$j,$k)=(0,0,0);
my @tm=();
my @dxy=();
my @fst=();

open(R,$in);
while(<R>)
 {chomp $_;
  @tm=split(/\t/,$_);
  if($tm[0] eq "Chr") {next;}
  my $snp_num=$tm[4];
  my $win_len=$tm[3]-$tm[2]+1;
  my @dxy_mat=@tm[$dxy_sta..$dxy_end];
  my @fst_mat=@tm[$fst_sta..$fst_end];
  if($snp_num>10)
    { open(WF,">infile");	
      printf WF "%5d\n",$sn;
      $k=0;
      for $i(0..$sn-2) { for $j($i+1..$sn-1) {$dxy[$i][$j]=$dxy_mat[$k]; $dxy[$j][$i]=$dxy[$i][$j];$k++;}} 
      for $i(0..$sn-1) 
            {  printf WF "%-10s",$sample[$i];
               for $j(0..$sn-1) 
               {  if($i==$j){$dxy[$i][$j]=0.0;}
                  printf WF "  %8.6f",$dxy[$i][$j];
                }
               print WF "\n";
             } 
      close(WF);
      system ("./neighbor");
      my $dtree="";
      open(RDT, "outtree");
      while(<RDT>){chomp $_; $dtree .= $_;} 
      close(RDT);
      printf WDT "%s\t%s\t%s\t%s\t%s\n",$tm[0],$tm[1],$tm[2],$tm[3],$dtree;
      unlink("infile");
      unlink("outfile");
      unlink("outtree");
      
      open(WF,">infile");	
      printf WF "%5d\n",$sn;
      $k=0;
      for $i(0..$sn-2) { for $j($i+1..$sn-1) {$fst[$i][$j]=$fst_mat[$k]; $fst[$j][$i]=$fst[$i][$j];$k++;}} 
      for $i(0..$sn-1) 
            {  printf WF "%-10s",$sample[$i];
               for $j(0..$sn-1) 
               {  if($i==$j){$fst[$i][$j]= 0.0;}
                  printf WF "  %8.6f",$fst[$i][$j];
               }
               print WF "\n";
             } 
      close(WF);
      system ("./neighbor");
      my $ftree="";
      open(RFT, "outtree");
      while(<RFT>) {chomp $_; $ftree .= $_;}
      close(RFT);
      printf WFT "%s\t%s\t%s\t%s\t%s\n",$tm[0],$tm[1],$tm[2],$tm[3],$ftree;
      unlink("infile");
      unlink("outfile");
      unlink("outtree");
     }
}
close(WDT);
close(WFT);
  
