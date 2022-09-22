#!/usr/bin/perl -w
my $snp_fre  = shift;  # allele frequency of SNPs 
my $snp_pos  = shift;  # position index of SNPs
my $wnd_size = shift;  # sliding windows size(kb)
my $wnd_step = shift;  # sliding step(kb)
my $out      = shift;  # output file

my @cr=qw(chr01 chr02 chr03 chr04 chr05 chr06 chr07 chr08 chr09 chr10 chr11 chr12);  

my @tm=();
my @tn=();
my ($i, $j);

my %pos_idx=();
open(R,$snp_pos) || die "reading $gene_info failed ..\n";
while( <R> )
     { chomp $_;
       @tm=split(/\t/, $_);
       $pos_idx{$tm[0]} = $tm[1];
     }
close(R);     

open(R, $snp_fre) || die "reading $gene_info failed ..\n";
my $headline = <R>;
chomp $headline; 
@tm=split(/\t/, $headline);
my @pop_name=qw(cultivated wild);
my %pop_fr=();
my %pop_pos = ();
for $i(5..scalar(@tm)-1) {$pop_pos{$tm[$i]} = $i;}
$i=0;

open(W ,">".$out) || die " writing $out failed ..\n";
printf W "%s\t%s\t%s\t%s\t%s", 'chr','segment','start_pos','end_pos','snp_num';
foreach(@pop_name) { printf W "\t%s-pi(1kb)", $_;
                     printf W "\t%s-theta-W(1kb)", $_; 
                     printf W "\t%s-Tajima-D\/sample", $_;
                    }
printf W "\t%s-%s-dxy(1kb)", $pop_name[0], $pop_name[1];
printf W "\t%s-%s-fst", $pop_name[0], $pop_name[1];
printf W "\t%s-%s-rod", $pop_name[0], $pop_name[1];
print W "\n";     
             

$wnd_size = $wnd_size*1000;  
$wnd_step = $wnd_step*1000;  
my ($start_pos,$end_pos);
foreach my $tcr(@cr)
       { my $flag=1;
       	 $start_pos=1;
       	 while($flag)
       	   { foreach(@pop_name){@{$pop_fr{$_}}=();}
             $end_pos = $start_pos + $wnd_size-1 ;
             my $len = $end_pos-$start_pos+1; 
             printf W "%s\t%s\t%s",$tcr,$start_pos,$end_pos;
             my $pos_marker=int($start_pos/10000);
             while(!exists($pos_idx{$tcr.'-'.$pos_marker})){$pos_marker--;}
             seek(R, $pos_idx{$tcr.'-'.$pos_marker},0);
             my $snm=0;
             while(<R>)
                { chomp $_;
     	          @tm=split(/\t/, $_);
                  if( $tm[0] eq $tcr && $start_pos <= $tm[1] && $tm[1] <= $end_pos )
                  { $snm++; foreach(@pop_name){ push(@{$pop_fr{$_}},$tm[$pop_pos{$_}]);}}         
                  if(($tm[0] eq $tcr && $tm[1] > $end_pos)) {last;}
                  if($tm[0] ne $tcr || eof(R)) {$flag=0; last}
                }
           printf W "\t%s", $snm;    
           my %k =();
           my %ss=();
           my %sn=();
           my %pi=();
           my %th=();
           my %tj=();
           foreach(@pop_name){ $k{$_}  = cal_k  (\@{$pop_fr{$_}});
                               $ss{$_} = cal_ss (\@{$pop_fr{$_}});
                               $sn{$_} = cal_sn (\@{$pop_fr{$_}});
                               $pi{$_} = 1000*$k{$_}/$len;
                               $th{$_} = 1000*cal_theta($ss{$_},$sn{$_},$len); 
                               $tj{$_} = cal_tajimaD($k{$_},$ss{$_},$sn{$_}); 
                              }
           my $dxy=0;
           my $fst=0;
           my $rod=0; 
           my $p1 = $pop_name[0];
           my $p2 = $pop_name[1];
           $dxy = cal_dxy(\@{$pop_fr{$p1}},\@{$pop_fr{$p2}},$len);
           $dxy = 1000*$dxy;
           if($dxy==0){ $fst=0;}
           else           { $fst=(2*$dxy-$pi{$p1}-$pi{$p2})/(2*$dxy); 
                            if($fst<0) {$fst=0;}
                          }
           if($pi{$p1}==0){$rod=999;}
           else           {$rod=$pi{$p2}/$pi{$p1};}
           foreach(@pop_name){ printf W "\t%10.5f", $pi{$_};
                               printf W "\t%10.5f", $th{$_}; 
                               printf W "\t%10.5f\/%d", $tj{$_},$sn{$_};
                             } 
           printf W "\t%10.5f", $dxy;
           printf W "\t%10.5f", $fst;
           printf W "\t%10.5f\n", $rod;
           $start_pos += $wnd_step; 
          }
      }                            
       

sub  cal_k  
       { my ($fr) = @_;
         my $k = 0;
         my @sfr = ();
         foreach (@$fr)
                 { @sfr=split(/:/,$_);
                   my $an=scalar(@sfr);
                   my ($i,$j);
                   my $mch=0;
                   for $i(0..$an-2) { for $j($i+1..$an-1) { $mch += $sfr[$i] * $sfr[$j]; } }
                   my $sun=0;
                   foreach(@sfr) { $sun +=$_;}
                   my $ek=0;
                   if($sun>1) {$ek=2*$mch/($sun*($sun-1));}
                   $k += $ek;
                  }
         return $k;
         }         

sub  cal_ss 
       { my ($fr) = @_;
         my $ss = 0;
         my @sfr = ();
         foreach (@$fr)
                 { @sfr=split(/:/,$_);
                   my $i=0;
                   foreach(@sfr){if($_ > 0) {$i++;}}
                   $ss += $i-1;
                 }
         return $ss;
        }  

sub  cal_sn
       { my ($fr) = @_;
         my $n=0;
         my $fn=0;
         my @sfr = ();
         foreach (@$fr)
                 { @sfr=split(/:/,$_);
                   my $i=0;
                   my $j=0;
                   foreach(@sfr){if($_ > 0) {$i++;}}
                   foreach(@sfr){$j += $_;}
                   if($i-1>0) {$n++; $fn += $j;}
                 }
         my $sn=0;
         if($n>0){$sn=int($fn/$n+0.5);}        
         return $sn;
        }  

sub cal_dxy 
      { my ($p1, $p2, $ln) = @_;
        my $dxy = 0;
        my @sfr1 = ();
	my @sfr2 = ();
	my ($i, $j);
	my $sn=scalar(@$p1);
	for $i(0..$sn-1)
	      { @sfr1=split(/:/,$p1->[$i]);               
                @sfr2=split(/:/,$p2->[$i]);
                my $mch=0;
                my $an=scalar(@sfr1);
                my $am1=0;
                my $am2=0;
                foreach(@sfr1){$am1 += $_;}
                foreach(@sfr2){$am2 += $_;}
                my $am=$am1*$am2;
                for $i(0..$an-1)
                      {for $j(0..$an-1)
                             { if($i==$j){next;}
                               else      {$mch += $sfr1[$i]*$sfr2[$j];}
                             }
                      }
                if($am>0){$dxy += $mch/$am;}
               }
        $dxy = $dxy/$ln;
        return $dxy;
       }             
                             


# calcute theta-W 
sub cal_theta 
      { my ($ss, $sn, $ln) = @_;   # $ss: number of segregating sites  $sn: sample size  $ln: sequence length 
        my $a1=0;
        my $i=0;
        my $theta=0;
        if($ss>0 && $sn>0)
          {    for $i(1..$sn-1) {$a1 += 1/$i;}
               $theta = $ss/$a1/$ln;
          }     
        return $theta;
       }

# calcute tajima-D 
sub cal_tajimaD
     {  my ($k, $ss, $sn) = @_;
        my ($a1, $b1, $c1, $e1) = (0, 0, 0, 0);
        my ($a2, $b2, $c2, $e2) = (0, 0, 0, 0);
        my ($i, $j);
        my $td=0;
        if($ss>1 && $sn>0)
          { for $i(1..$sn-1) {$a1 += 1/$i;}
            $b1 = ($sn+1)/(3*($sn-1));
            $c1 = $b1-1/$a1;
            $e1 = $c1/$a1;
            for $i(1..$sn-1) {$a2 += 1/($i*$i);}
            $b2 = 2*($sn*$sn+$sn+3)/(9*$sn*($sn-1));
            $c2 = $b2-($sn+2)/($a1*$sn)+$a2/($a1*$a1);
            $e2 = $c2/($a1*$a1+$a2);
            $td=($k-$ss/$a1)/sqrt($e1*$ss+$e2*$ss*($ss-1));
           }
        return $td;
      }         
        
        
       

     

       
                           
        
      
       
