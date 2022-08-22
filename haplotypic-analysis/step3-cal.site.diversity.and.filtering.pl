#!/usr/bin/perl -w

use strict;

my $gene   = shift;   
my $lbdr   = shift;
my $rbdr   = shift;
my $segdir = shift;
my $divdir = shift;


if(!opendir(DIV, $divdir)){ mkdir($divdir,oct(750)) || die "Establishing folder $divdir failed .. /n"; }
else                      { closedir(DIV); }

my %sam_grp=();
my @sam=();
my @tm=();
my ($i, $j);
open(R,"201110.sample.grp") || die "Couldn't open 200928.1493.samples.grp for reading ..\n";
while ( <R> )
      { chomp $_;
        @tm=split(/\t/, $_);
        push(@sam, $tm[0]);
        $sam_grp{$tm[0]}=$tm[1];
       }
 close(R);
my @corn=qw(IND TEJ TRJ);
my @cor=();
my @satn=qw(IND TEJ TRJ AUS ARO RAY);
my @sat=();
my @wildn=qw(RUF1 RUF2 NIV1 NIV2);
my @wild=();
my @pg=(\@wild, \@sat,\@cor);

foreach my $sm(@sam) {  foreach(@corn)   { if($_ eq $sam_grp{$sm})   {  push(@cor, $sm);}}
 	                foreach(@satn)   { if($_ eq $sam_grp{$sm})   {  push(@sat, $sm);}}
 	                foreach(@wildn)  { if($_ eq $sam_grp{$sm})   {  push(@wild, $sm);}}
 	             }


my @gen_mat = ();
open (R,$gene) || die "Couldn't open $gene for reading .. \n";
$i=0;
while(<R>)
     { chomp $_;
       $_ =~  s/ //g;
       @tm=split(/\t/,$_);
       push(@{$gen_mat[$i]},@tm);
       push(@{$gen_mat[$i]},'0');  
       
       $i++;
     }
close(R);

foreach my $gm(@gen_mat)
              { open(R,  $segdir."/". $gm->[0]. ".segment.vcf" ) || die "reading file fails ..\n";         
                my @sdt=<R>;
                close(R);
                my %sam_pos=();
                
                foreach (@sdt) {if ($_ =~ /^\#CHROM/ ){chomp $_; @tm=split(/\t/,  $_);last;}}
                for $i(0..8){shift(@tm);}
                $i=9; foreach(@tm) {$sam_pos{$_}=$i; $i++;}  
                
                my @fcsam=();
                foreach (@sdt)
                        { chomp $_;
                          if($_ =~ /^\#/) { next }
                          @tm=split(/\t/, $_);
                          if($tm[1]==$gm->[4])
                            {foreach(@cor)
                                    {if($tm[$sam_pos{$_}] =~ /^$gm->[5]/)
                                       {push(@fcsam, $_);}
                                    }
                             last;
                             }
                         }  
                
                open(W,  ">".$divdir."/". $gm->[0]. ".segment.pi" ) || die "writing file fails ..\n"; 
                print W "chr\tpos\tref\talt\tqual\twild_sn\twild_an\twild_fr\twild_pi\twild_he\tsat_sn\tsat_an\tsat_fr\tsat_pi\tsat_he\tcor_sn\tcor_an\tcor_fr\tcor_pi\tcor_he\tfcor_sn\tfcor_an\tfcor_fr\tfcor_pi\tfcor_he\tfst\n";
                foreach (@sdt)
                        { chomp $_;
                          if($_ =~ /^\#/) { next }
                          @tm=split(/\t/, $_);
                          my @allele=split(/,/,$tm[4]);
                          my $aln=scalar(@allele);
                          printf W "%s\t%s\t%s\t%s\t%s", $tm[0], $tm[1], $tm[3], $tm[4],$tm[5];
                          
                          my @sg=();
                          my  ($sn, $ps, $cn, $pi, $he);
                          my @fr=();
                          foreach(@wild){push(@sg, $tm[$sam_pos{$_}]);} 
                          $sn=scalar(@sg);
                          @fr=cal_fr($aln, \@sg);
                          $cn=0;
                          foreach(@fr) { $cn += $_;}
                          $cn=$cn/2;
                          $ps=$fr[0];
                          for $i(1..$aln)  { $ps .= "|".$fr[$i]; }  
                          $pi=cal_pi(\@fr);
                          $he=cal_he(\@fr);
                          my @wfr=@fr;
                          printf W "\t%s\t%s\t%s\t%8.5f\t%8.5f", $sn, $cn, $ps, $pi, $he;
                          
                          @sg=();
                          foreach(@sat) { push(@sg, $tm[$sam_pos{$_}]);} 
                          $sn=scalar(@sg);
                          @fr=cal_fr($aln, \@sg);
                          $cn=0;
                          foreach(@fr) {$cn += $_;}
                          $cn=$cn/2;
                          $ps=$fr[0];
                          for $i(1..$aln)  { $ps .= "|".$fr[$i];}          
                          $pi=cal_pi(\@fr);
                          $he=cal_he(\@fr);
                          printf W "\t%s\t%s\t%s\t%8.5f\t%8.5f", $sn, $cn, $ps, $pi, $he;
                         
                          
                          @sg=();
                          foreach(@cor) {  push(@sg, $tm[$sam_pos{$_}]);} 
                          $sn=scalar(@sg);
                          @fr=cal_fr($aln, \@sg);
                          $cn=0;
                          foreach(@fr) {$cn += $_;}
                          $cn=$cn/2;
                          $ps=$fr[0];
                          for $i(1..$aln)  { $ps .= "|".$fr[$i];}          
                          $pi=cal_pi(\@fr);
                          $he=cal_he(\@fr);
                          printf W "\t%s\t%s\t%s\t%8.5f\t%8.5f", $sn, $cn, $ps, $pi, $he;
                                                     
                          @sg=();
                          foreach(@fcsam) {  push(@sg, $tm[$sam_pos{$_}]);} 
                          $sn=scalar(@sg);
                          @fr=cal_fr($aln, \@sg);
                          $cn=0;
                          foreach(@fr) {$cn += $_;}
                          $cn=$cn/2;
                          $ps=$fr[0];
                          for $i(1..$aln)  { $ps .= "|".$fr[$i];}          
                          $pi=cal_pi(\@fr);
                          $he=cal_he(\@fr);
                          printf W "\t%s\t%s\t%s\t%8.5f\t%8.5f", $sn, $cn, $ps, $pi, $he;
                          my $fst=cal_fst(\@fr,\@wfr);
                          printf W "\t%8.5f\n", $fst; 
                         }
                close(W);

                open(R,  $divdir."/". $gm->[0]. ".segment.pi" ) || die "reading file fails ..\n"; 
                my @vdt=<R>;
                close(R);         
                open(W, ">".$divdir."/". $gm->[0]. ".filter.segment.pi" ) || die "writing file fails ..\n";
                foreach(@vdt) { chomp $_; }
                my $fl1=shift(@vdt);
                printf W "%s\t%s\t%s\n", $fl1, "MAF", "Missing_rate";
                foreach(@vdt)
                       { @tm=split(/\t/, $_);
                         if($tm[1] != $gm->[4] && ($tm[6]+$tm[11])/($tm[5]+$tm[10])<=0.9){next;}  
                         if(cal_maf($tm[7],$tm[12])<=0.05) {next;}
                         if($tm[3] =~ /\,/) { next;}
                         printf W "%s\t%8.5f\t%8.5f\n", $_, cal_maf($tm[7],$tm[12]), ($tm[6]+$tm[11])/($tm[5]+$tm[10]);
                        }
                close(W);
                
                
                
                open(R, $divdir."/". $gm->[0]. ".filter.segment.pi") || die "reading file fails ..\n"; 
                my @wdt=<R>;
                close(R);  

                open(W, ">".$divdir."/". $gm->[0]. ".1k.pi.filter" ) || die "writing file fails ..\n";
                my $fl2=shift(@wdt);
                print W $fl2;
                my $lbs=$gm->[2]-$lbdr; if($lbs<0){$lbs=0;}
                my $rbs=$gm->[3]+$rbdr; 
                foreach(@wdt)
                       { @tm=split(/\t/, $_);
                         if($lbs<=$tm[1] && $tm[1]<=$rbs){print W $_;}
                         if($tm[1]>$rbs){last;}
                       }
                close(W);       
               }

sub   cal_fr  {   my ($an, $seg)=@_;
                         my %ale=();
                         my $i=0;
                         for $i(0..$an) {$ale{$i}=0;}
                         foreach (@$seg)
                                      {  if( $_ =~ /^\./) {next;} 
                                         my @tk=split(/:/,$_);
                                         my @ta=split(/\//, $tk[0]);
                                         $ale{$ta[0]}++; 
                                         $ale{$ta[1]}++;
                                      }
                         my @sfr=();
                         for  $i(0..$an) { push (@sfr, $ale{$i}); }
                         return @sfr;
                       }

sub  cal_he     {  my ($sfr) = @_;
                       my $an=scalar(@$sfr);
                       my $atn=0;
                       my @afr=();
                       my $af=0;
                       my $he=0;
                       foreach(@$sfr){$atn += $_;}
                       if($atn>0) { foreach(@$sfr){$af=$_/$atn; push(@afr, $af);}
                                    my $sq=0;
                                    foreach(@afr){ $sq += $_*$_;}
                                    $he=1-$sq;
                                  }
                       return $he;
                   }



sub  cal_pi  {  my ($sfr) = @_;
                       my $an=scalar(@$sfr);
                       my ($i,$j);
                       my $mch=0;
                       for $i(0..$an-2) { for $j($i+1..$an-1) {$mch += $sfr->[$i] * $sfr->[$j]; }}
                       my $sun=0;
                       foreach(@$sfr) { $sun +=$_;}
                       my $pi=0;
                       if($sun>1) { $pi=2*$mch/($sun*($sun-1));}
                       return $pi;
                      }

sub cal_fst  { my ($sfr1, $sfr2)=@_;
                      my $an=scalar(@$sfr1); 
                      my ($i,$j);
                       my $mch=0;
                       for $i(0..$an-1) {for $j(0..$an-1) { if($i!=$j) { $mch +=$sfr1->[$i] *$sfr2->[$j]; }}}
                       my $sun1=0;  foreach(@$sfr1) { $sun1 +=$_;}
                       my $sun2=0;  foreach(@$sfr2) { $sun2 +=$_;}
                       my $dxy=0;
                       if( $sun1>0  && $sun2>0) { $dxy=$mch/($sun1*$sun2); }
                       my $pi1=cal_pi($sfr1);
                       my $pi2=cal_pi($sfr2);
                       my $fst=0;
                       if ($dxy>0) { $fst=($dxy-0.5*($pi1+$pi2))/$dxy;}
                       if($fst<0) { $fst=0;}
                       return $fst;
                     }


sub cal_maf { my ($f1, $f2)=@_;
              my @fr1=split(/\|/, $f1);
     	      my @fr2=split(/\|/, $f2);
     	      my $fn=scalar(@fr1);
     	      my @fr=();
     	      my $i=0;
              my $mxf=0;
              my $smf=0;
              for $i(0..$fn-1) { $fr[$i] =$fr1[$i]+$fr2[$i]; 
                                 $smf += $fr[$i];
                                 if($fr[$i]> $mxf){$mxf=$fr[$i];}
                               }
              my $maf=0;
              if($smf>0) {$maf=($smf-$mxf)/$smf;}
              return $maf;
             }                  
              
              
     	              
	
                                                                                             	             
