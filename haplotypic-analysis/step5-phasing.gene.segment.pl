#The java program beagle.12Jul19.0df.jar used in this perl script is from the software BEAGLE      

#!/usr/bin/perl -w

use strict;

my $genlst=shift;
my $segdir=shift;
my $phsdir=shift;

if(!opendir(PHS, $phsdir)){ mkdir($phsdir,oct(750)) || die "Establishing folder $phsdir failed .. /n"; }
else                      { closedir(PHS); }


my @gene = ();
my @tm = ();
open(R, $genlst) || die "can not open $genlst ../n";
while ( <R> )  {@tm=split(/\t/, $_);   push(@gene,$tm[0]); }
close(R);

foreach my $gn(@gene)  
        {  my $fn= $gn . ".segment.vcf.gz";                                     
	   my $gtf=$fn;
           $gtf  =~ s/.vcf.gz//;
           $gtf .= '.gt';
           system ("java -jar beagle.12Jul19.0df.jar gt=$segdir/$fn out=$phsdir/$gtf");
           my $fpn = $segdir."/".$gn.".1k.pi.filter";
           open(R,$fpn) || die "reading failed..\n";
           my %rst=();
           while(<R>) {@tm=split(/\t/, $_); $rst{$tm[1]}=1; }  
           close(R);     
           my $ftn=$phsdir."/".$gn.".segment.gt.vcf.gz";
           my $wtn=$phsdir."/".$gn.".hap.vcf";
           open(R,"gzip -dc $ftn|") || die "reading failed..\n";
           open(W,">".$wtn) || die "writing failed..\n";
           while(<R>) 
                { if ($_ =~ /^\#/) {print W $_; next;}
                  @tm = split(/\t/, $_);
                  if(exists($rst{$tm[1]})){print W $_;}
                }   
           close(R);
           close(W);   
         }
                 
                   
