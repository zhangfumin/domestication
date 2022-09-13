#!/usr/bin/perl -w

use strict;

my $gene_info  = shift;  # the file containing gene infomation. Data format: gene_name        chr     start_pos       tend_pos
                         # ie:   OS100000       chr01   100345  101765
my $upstream   = shift;  # upstream, bp 
my $downstream = shift;  # downstream, bp
my $bootstrip  = shift;  # replication number

srand(time());

my ($i,$j,$b);
my @tm=();

open(R,$gene_info) || die "Failed to load $gene_info, please check the file $gene_info."  ; 
my @g_info=<R>; 
close(R);

foreach (@g_info) 
        { chomp $_; 
          my ($gn, $cr, $stp, $edp) = split(/\t/, $_);
          my $fre = "group_".$cr.".fr";   #frequency file of each chromosome calculated by the program cal.pop.fr.pl
          my $idx = $fre.".idx";           #frequency file index built by the program build.fr.idx.pl
          my %pos_idx=();
          open (R, $idx) || die "reading failed..\n";
          while (<R>){ chomp $_; @tm=split(/\t/,$_); $pos_idx{$tm[0]} = $tm[1]; }
          close(R); 
          $stp = $stp - $upstream;  if($stp<0) { $stp = 0;}
          $edp = $edp + $downstream;
          my $pos_marker=int($stp/10000);
          while(!exists($pos_idx{$cr.'-'.$pos_marker})){$pos_marker--;}
          open ( R, $fre) || die "reading failed..\n";
          open ( W, ">".$gn.".".$bootstrip.".bt") || die "writing failed..\n";
          my $headline=<R>;
          chomp $headline;
          @tm=split(/\t/, $headline);
          my @subpop_name=@tm;
          for $i(0..3) {shift(@subpop_name);}
          printf W "%s\t%s\t%s\t%s\t%s", "Gene","Chr","Start_pos","End_pos","Snp_num";
          foreach (@subpop_name){printf W "\tPi_%s", $_;} 
          my $sn=scalar(@subpop_name);   $sn--;
          for $i(0..$sn-1) {for $j($i+1..$sn) {printf W "\tDxy_%s-%s", $subpop_name[$i],$subpop_name[$j];}}
          for $i(0..$sn-1) {for $j($i+1..$sn) {printf W "\tFst_%s-%s", $subpop_name[$i],$subpop_name[$j];}}
          print W "\n";
          seek(R, $pos_idx{$cr.'-'.$pos_marker},0);
          my @win=();
          my $gr=1; 
          while( $gr )
               { my $line=<R>;
                 chomp $line;
                 @tm=split(/\t/, $line);    
                 if( $stp <= $tm[1] && $tm[1]<=$edp ) {push(@win, $line);}
                 if($tm[1]>$edp) {$gr=0;}
               }
          close(R);
          my $wn=scalar(@win);
          for $b(0..$bootstrip-1 )
                { my @bw=();
                  for $j(0..$wn-1){push(@bw, int(rand($wn)));}
                  my @pi  = ();
                  my @dxy = ();
                  my @fst = ();
                  for $i(0..$sn) { $pi[$i]=0;}
                  for $i(0..$sn-1) { for $j($i+1..$sn) {$dxy[$i][$j]=0;}}
                  for $i(0..$sn-1) { for $j($i+1..$sn) {$fst[$i][$j]=0;}}
                  foreach (@bw) { @tm=split(/\t/, $win[$_]);
                                  for $i(0..3)   { shift(@tm);}
                                  for $i(0..$sn) { $pi[$i] += cal_pi($tm[$i]); }
                                  for $i(0..$sn-1) { for $j($i+1..$sn) { $dxy[$i][$j] += cal_dxy($tm[$i],$tm[$j]);}}
                                }
                   for $i(0..$sn-1){ for $j($i+1..$sn) 
                                           { if($dxy[$i][$j]!=0){$fst[$i][$j]=1-($pi[$i]+$pi[$j])/(2*$dxy[$i][$j]);}
                                             if($fst[$i][$j]<0) {$fst[$i][$j]=0;}
                                            }
                                    }
                   my $snp_num=scalar(@bw);
                   my $start_pos = $stp;
                   my $end_pos   = $edp;
                   my $win_len=$end_pos-$start_pos+1; 
                   printf W "%s\t%s\t%s\t%s\t%s",$gn,$cr,$start_pos,$end_pos,$snp_num;
                   for $i(0..$sn)   { printf W "\t%8.3f",1000*$pi[$i]/$win_len; }
                   for $i(0..$sn-1) { for $j($i+1..$sn) {printf W "\t%8.3f",1000*$dxy[$i][$j]/$win_len;}}
                   for $i(0..$sn-1) { for $j($i+1..$sn) {printf W "\t%8.6f",$fst[$i][$j];}}                  
                   print W "\n";
                  }
               close(W);
             }   
                                   
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
         
