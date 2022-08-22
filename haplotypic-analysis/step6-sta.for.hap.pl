#!/usr/bin/perl -w

use strict;

my $genlst=shift;
my $phsdir=shift;
my $stadir=shift;

if(!opendir(STA, $stadir)){ mkdir($stadir,oct(750)) || die "Establishing folder $stadir failed .. /n"; }
else                      { closedir(STA); }

open(R, $genlst) || die "can not open $genlst ../n";
my @gene=();
my @tm = ();
my ($i,$j);
while ( <R> ) {@tm=split(/\t/, $_); push(@gene,$tm[0]);}
close(R);

my @gpn = qw(IND TEJ TRJ RUF1 RUF2 NIV1 NIV2);

open(R,"201110.sample.grp") || die "Couldn't open 200928.1493.samples.grp for reading ..\n";  
my %grp=();
while ( <R> ){chomp $_;  @tm=split(/\t/, $_); $grp{$tm[0]}=$tm[1];}
close(R);      

foreach my $gnm(@gene)
         { my $fn_gene =  $phsdir.'/'.$gnm.'.hap.vcf';
           my $fn_hap  =  $stadir.'/'.$gnm.'.hap.sta';
           open(R, $fn_gene)    || die "Cannot open $fn_gene for reading ..";
           open(W, ">".$fn_hap) || die "Cannot write $fn_hap .."; 
           my $headline=0;
           my @tsam=();          
           my @sam=();
           my @pdt=();
           while(<R>) 
               { if($_ =~ /^\#CHROM/) {$headline = $_;}
                 if($_ =~ /^\#/)      {next;}
                 chomp $_;
                 push(@pdt, $_);
               } 
          close(R);  
          my @pgp=();
          for $i(0..scalar(@gpn)-1){@{$pgp[$i]}=();}
          my %sam_pos = ();
          chomp $headline;
          @tm=split(/\t/, $headline); 
          for $i(0..8) { shift(@tm);}
          $i=9; foreach (@tm){push(@tsam, $_); $sam_pos{$_}=$i; $i++;}
          foreach my $tsm(@tsam)
                    { my $fg=0;
                      foreach(@gpn) {if($grp{$tsm} eq $_) {$fg=1; last;}}
                      if($fg) {push(@sam, $tsm);}
                     } 
          foreach $j(@sam) {for $i(0..scalar(@gpn)-1){if($grp{$j} eq $gpn[$i]){push(@{$pgp[$i]},$j);}}} 
          my @snp_pos=();
          my %hap_phs=();
          my %het_num=();
          my @sit_het=();
          my @tn=();
          my $snp=scalar(@pdt);
          foreach(@sam)    { $hap_phs{$_.'.1'}=""; $hap_phs{$_.'.2'}=""; $het_num{$_}=0; }
          for $i(0..$snp-1){ $sit_het[$i]=0;}
          $i=0;     
          my @var=();
          foreach (@pdt) { @tm=split(/\t/, $_);
                           push(@var, sprintf("%s\t %s\t%s", $tm[1],$tm[3],$tm[4]));
                           push(@snp_pos, $tm[1]);
                           foreach(@sam)
                                  { @tn=split(/\|/, $tm[$sam_pos{$_}]); 
                                    $hap_phs{$_.'.1'} .= $tn[0]; 
                                    $hap_phs{$_.'.2'} .= $tn[1];
                                    if($tn[0] ne $tn[1]) {$het_num{$_}++; $sit_het[$i]++;} 
                                  }
                           $i++;
                          }
          my @hap_dat=();
          foreach(@sam) { push(@hap_dat, $hap_phs{$_.'.1'}); push(@hap_dat, $hap_phs{$_.'.2'});}
          my @hap=();
          push(@hap, $hap_dat[0]);
          my $dat=0;
          my $flag=0;
          foreach $dat( @hap_dat )
                      { $flag=1;
                        foreach (@hap) {if($dat eq $_) {$flag=0; last;}}
                        if($flag){push(@hap,$dat);}
                       }
         
          my $hap_total=scalar(@hap_dat);
          my %hap_num=();
          foreach(@hap)     {$hap_num{$_}=0;}
          foreach(@hap_dat) {$hap_num{$_}++;}
          my @nhp=();
          foreach(@hap){if($hap_num{$_}/$hap_total>0.005){push(@nhp, $_);}}
          
          my %nhp_no=();
          $i=1;
          foreach(@nhp) {$nhp_no{$_} = 'H'.$i; $i++;} 
          
          printf W "Variation of  sites\n";
          $i=1;
          foreach (@var)  { printf W  "%s\t%s\n",$i, $_; $i++;}
          print W "\n";              
          printf W "Statistics for %s haplotype\n", $gnm;
          printf W "A total of %d %s haplotypes in rice\n", scalar(@nhp), $gnm;    
          $i=1; foreach(@nhp){ printf W "H%s\t%s\n", $i, $_; $i++;}
          print W "\n"; 
          
          printf W "The distribution of %s haplotype in rice groups\n", $gnm;
          print W "Haplotype\tAll";
          foreach(@gpn){ printf W "\t%s", $_;}
          print W "\n"; 
          $i=1;
          foreach my $np(@nhp)
                     { printf W "H%s", $i;
                       my $htn=0;
                       foreach(@sam) { if(exists($nhp_no{$hap_phs{$_.'.1'}}) && exists($nhp_no{$hap_phs{$_.'.2'}}))
                                         { if($hap_phs{$_.'.1'} eq $np) {$htn++;}
                                           if($hap_phs{$_.'.2'} eq $np) {$htn++;}
                                         }
                                      }
                       printf W "\t%s", $htn;                  
                                           
                       for $i(0..scalar(@gpn)-1)
                          { $htn=0;
                            foreach(@{$pgp[$i]})
                                    { if(exists($nhp_no{$hap_phs{$_.'.1'}}) && exists($nhp_no{$hap_phs{$_.'.2'}}))
                                        { if($hap_phs{$_.'.1'} eq $np) {$htn++;}
                                          if($hap_phs{$_.'.2'} eq $np) {$htn++;}
                                        }  
                                    }
                             printf W "\t%s", $htn;
                           }
                          print W "\n";
                        $i++;
                     }
          print W "\n";
          printf W "The hyplotype of each sample in each group\n";
          
          
          
          for $i(0..scalar(@gpn)-1)
                {  printf W "%s\t%s\n", $gpn[$i], "haplotype";
                   foreach(@{$pgp[$i]}) { if(exists($nhp_no{$hap_phs{$_.'.1'}}) && exists($nhp_no{$hap_phs{$_.'.2'}}))
                   	                  { printf W "%s_1\t%s\n", $_,$nhp_no{$hap_phs{$_.'.1'}};
                                            printf W "%s_2\t%s\n", $_,$nhp_no{$hap_phs{$_.'.2'}};
                                          }
                                        } 
                }
          close(W);
        }                                       
                   
