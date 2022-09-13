#!/usr/bin/perl -w

use strict;

my $in = shift;  # A file containing SNP allele frequency of populations calculated by the program cal.pop.fr.pl 
my $gene_info = shift; # A file containing windows or genes infomation. Every line includes 4 fields, window or gene name,chromosome, start position and end position, which were separated by tab characters. For example, OS100000	chr01	100345	101765
my $out = shift;  # An output file for storing results

my ($i,$j);
my @tm=();
my @tn=();

open(R,$in) || die "Failed to load $in, please check the file $in."  ; 
my @dt=<R>; 
close(R);
foreach (@dt) { chomp $_; }

open(R,$gene_info) || die "Failed to load $gene_info, please check the file $gene_info."  ; 
my @gn=<R>; 
close(R);
foreach (@gn) { chomp $_; }

my @gnw=();
my $gsn=scalar(@gn);

for $i(0..$gsn-1){ @{$gnw[$i]}=(); }

@tm=split(/\t/, $dt[0]);
my @subpop_name=@tm;
for $i(0..3) {shift(@subpop_name);}

open(W,">".$out);
printf W "%s\t%s\t%s\t%s\t%s", "Gene","Chr","Start_pos","End_pos","Snp_num";
foreach (@subpop_name){printf W "\tPi_%s", $_;} 
my $sn=scalar(@subpop_name);   $sn--;
for $i(0..$sn-1) {for $j($i+1..$sn) {printf W "\tDxy_%s-%s", $subpop_name[$i],$subpop_name[$j];}}
for $i(0..$sn-1) {for $j($i+1..$sn) {printf W "\tFst_%s-%s", $subpop_name[$i],$subpop_name[$j];}}
print W "\n";
shift(@dt);

my @win=();

foreach $j(@dt)
    { @tm=split(/\t/, $j);    
      for $i(0..$gsn-1){ @tn=split(/\t/,$gn[$i]);
                         push(@{$gnw[$i]},$j) if($tm[0] eq $tn[1] && $tn[2]<=$tm[1] && $tm[1]<=$tn[3]);                         ;     
                       }
    }                     
  for $i(0..$gsn-1)
   { @tn=split(/\t/,$gn[$i]);
     @win=@{$gnw[$i]};
     my @pi  = ();
     my @dxy = ();
     my @fst = ();
     for $i(0..$sn) { $pi[$i]=0;}
     for $i(0..$sn-1) { for $j($i+1..$sn) {$dxy[$i][$j]=0;}}
     for $i(0..$sn-1) { for $j($i+1..$sn) {$fst[$i][$j]=0;}}
     foreach (@win) { @tm=split(/\t/,$_);
                      for $i(0..3)     { shift(@tm);}
                      for $i(0..$sn)   { $pi[$i] += cal_pi($tm[$i]); }
                      for $i(0..$sn-1) { for $j($i+1..$sn) { $dxy[$i][$j] += cal_dxy($tm[$i],$tm[$j]);}}
                    }
    for $i(0..$sn-1){ for $j($i+1..$sn) {if($dxy[$i][$j]!=0){$fst[$i][$j]=1-($pi[$i]+$pi[$j])/(2*$dxy[$i][$j]);}
                                         if($fst[$i][$j]<0) {$fst[$i][$j]=0;}
                                       }
                    }
    my $snp_num=scalar(@win);
    my $start_pos = $tn[2];
    my $end_pos   = $tn[3];
    my $win_len=$end_pos-$start_pos+1; 

    printf W "%s\t%s\t%s\t%s\t%s",$tn[0],$tn[1],$start_pos,$end_pos,$snp_num;
    for $i(0..$sn)   { printf W "\t%8.3f",1000*$pi[$i]/$win_len; }
    for $i(0..$sn-1) { for $j($i+1..$sn) {printf W "\t%8.3f",1000*$dxy[$i][$j]/$win_len;}}
    for $i(0..$sn-1) { for $j($i+1..$sn) {printf W "\t%8.6f",$fst[$i][$j];}}                  
    print W "\n";
   }
    
close(W);                                        
#======================================================================================

sub cal_pi
    { my ($mp)=@_;
      my @tn=split(/:/,$mp);
      my $total=0;
      my $mismatch=0;
      my $pi=0;
      my ($c,$k,$an);
      $an=scalar(@tn); $an--;
      foreach(@tn){$total += $_;}
      for $c(0..$an-1) {for $k($c+1..$an) {$mismatch += $tn[$c]*$tn[$k];}}
      if($total > 1){ $pi=2*$mismatch/($total*($total-1)); }
      return $pi;
    }

sub cal_dxy
   { my($mp,$np)=@_;
     my @tn1=split(/:/,$mp);
     my @tn2=split(/:/,$np);
     my $total1=0;
     my $total2=0;
     my $dxy=0;
     foreach(@tn1){$total1 += $_;}
     foreach(@tn2){$total2 += $_;}
     my $total = $total1*$total2;
     if($total > 0)
       {  my $match=0;
          my ($c,$an);
          $an = scalar(@tn1);  $an--;   
          for $c(0..$an){ $match += $tn1[$c]*$tn2[$c];}
          my $mismatch=$total-$match;
          $dxy=$mismatch/$total;
      }
     return $dxy;
    }          
         
