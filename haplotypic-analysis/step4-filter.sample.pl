#!/usr/bin/perl -w

use strict;

my $gene   = shift;   
my $segdir = shift;
my $divdir = shift;
my $nsgdir = shift;

if(!opendir(NSG, $nsgdir)){ mkdir($nsgdir,oct(750)) || die "Establishing folder $nsgdir failed .. /n"; }
else                      { closedir(NSG); }

my @gen_mat = ();
open (R,$gene) || die "Couldn't open $gene for reading .. \n";
my @tm=();
my $i=0;
while(<R>)
     { chomp $_;
       $_ =~  s/ //g;
       @tm=split(/\t/,$_);
       push(@gen_mat,$tm[0]);
     }
close(R);

foreach (@gen_mat)
           { my $fseg  = $segdir.'/' .$_. '.segment.vcf';  
             my $fsit  = $divdir.'/' .$_. '.1k.pi.filter';
             my $rsit  = $divdir.'/' .$_. '.filter.segment.pi';
             my $nseg  = $nsgdir.'/' .$_. '.segment.vcf';    
             my $smis  = $nsgdir.'/' .$_. '.sample.missing'; 
             my $ssit  = $nsgdir.'/' .$_. '.1k.pi.filter';       

             my %retain=();
             my $sitn=0;
             open(R,$fsit) || die "Cannot open $fsit for reading ..";
             open(W,">".$ssit) || die "Cannot open $ssit for writing ..";
             <R>;
             while(<R>)
                  { print W $_;
                    @tm=split(/\t/,$_); $retain{$tm[1]} = 1;  $sitn++;
                  } 
             close(R);          
             close(W);
             
             my %rst=();
             open(R,$rsit) || die "Cannot open $fsit for reading ..\n";
             my $fl=1;
             while(<R>) { if($fl){$fl=0; next;}
             	          @tm=split(/\t/,$_); 
                          $rst{$tm[1]}=1;
                        }
             close(R);
             close(W);
             my %mis=();
             my @mdt=();
             my @msm=();
             my %msp=();
             my %msd=();
             open ( R, $fseg )    || die "Cannot open $fseg for reading ..";
             while(<R>) 
                  {  if($_ =~ /^\#CHROM/) { chomp $_; @tm=split(/\t/, $_); next;}
                     if($_ =~ /^\#/) {next;}
                     chomp $_;
                     push(@mdt, $_);
                  } 
             close(R);  
             for $i(0..8)  { shift(@tm); }
             $i=9; foreach (@tm){ push(@msm,$_);  $msp{$_}=$i; $i++;}
             foreach(@msm) {$msd {$_}=0;}
             foreach (@mdt)
                     { @tm=split(/\t/, $_);
                       if ($retain{$tm[1]}) {foreach(@msm)  { if( $tm[$msp{$_}] =~ /^\./) { $msd{$_}++;} } }
                     }
             my @rsm=();
             open(W,">".$smis ) || die " writing file failed.. \n";
             foreach(@msm) { printf W "%s\t%s\t%s\t%8.3f\n", $_, $msd{$_}, $sitn, $msd{$_}/$sitn;
                             if($msd{$_}/$sitn<=0.2) { push(@rsm, $_);}
                           }     
             close(W);
             
             open ( R, $fseg )    || die "Cannot open $fseg for reading ..";
             open ( W,">".$nseg)    || die " writing file failed.. \n";
             my $hf=0;
             while(<R>)
                  { if($_ =~ /^\#\#/) {print W $_; next;}
                    chomp $_;
                    @tm=split(/\t/, $_);
                    if($_ =~ /^\#CHROM/) {$hf=1;}
                    else                 {$hf=0;}
                    if($hf||exists($rst{$tm[1]})) 
                       {printf W "%s",$tm[0];
                        for $i(1..8){printf W "\t%s",$tm[$i];}
                        foreach(@rsm){printf W "\t%s",$tm[$msp{$_}];}
                        print W "\n";
                       } 
                  }
            close(R);
            close(W);
            system("gzip $nseg");
           }
           
                 
                            
             
                                                                                                 	             
