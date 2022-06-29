#Demography analysis for Oryza rufipogon (outcrossing) samples

    #Generating whole-genome diploid consensus sequence
    samtools mpileup -uf <reference> <sample>.bam | bcftools view -c - | vcfutils.pl vcf2fq -d <min_depth> -D <max_depth> | gzip > <sample>.diploid.fq.gz

    #Transforming the consensus sequence into the input file of psmc
    fq2psmcfa -q20 <sample>.diploid.fq.gz > <sample>.diploid.psmcfa

    #Inferring the population size history using 'psmc' (PSMC model by Heng Li) 
    psmc -N20 -t30 -r5 -p "4+30*2+4+6+10" -o <sample>.diploid.psmc <sample>.diploid.psmcfa


#Demography analysis for Oryza nivara and cultivate rice (selfing) samples

    #Generating VCF and Mask files from individual bam files
    samtools mpileup -q 20 -Q 20 -C 50 -u -r <chr> -f <reference> <sample1>.bam | bcftools view -cgI - | bamCaller.py <mean_depth> <sample1>_<chr>_mask.bed.gz | gzip -c > <sample1>_<chr>.vcf.gz
    samtools mpileup -q 20 -Q 20 -C 50 -u -r <chr> -f <reference> <sample2>.bam | bcftools view -cgI - | bamCaller.py <mean_depth> <sample2>_<chr>_mask.bed.gz | gzip -c > <sample2>_<chr>.vcf.gz

    #Combining two individuals into one input file
    generate_multihetsep.py --mask=<sample1>_<chr>_mask.bed.gz --mask=<sample2>_<chr>_mask.bed.gz <sample1>_<chr>.vcf.gz <sample2>_<chr>.vcf.gz > <sample1>-<sample2>.<chr>.msmc2.input

    #For each individual, randomly chosing one allele at occasional heterozygous sites to create pseudodiploid genomes (see Supplementary Materials 1.3)
    awk 'BEGIN{
        srand(19050801)
    }
    {
    if($4!~/,/){
        print $0
    } else {
        split(substr($4,1,4),allel,"")
        
        if(allel[1]!=allel[2]){
            a=allel[((rand()>=0.5)+1)]
        } else {
            a=allel[1]
        }

        if(allel[3]!=allel[4]){
            b=allel[((rand()>=0.5)+3)]
        } else {
            b=allel[3]
        }

        print $1"\t"$2"\t"$3"\t"a""a""b""b
    }
    }' <sample1>-<sample2>.<chr>.msmc2.input > <sample1>-<sample2>.<chr>.msmc2.input.rand_allel

    #Inferring the population size history using 'msmc2' (PSMC' model, a new implementation of PSMC model by Stephan Schiffels)
    msmc2 -t 1 -p 1*2+25*1+1*2+1*3 -o <sample1>-<sample2>.msmc2 -I 0-2 <sample1>-<sample2>.<chr*>.msmc2.input.rand_allel

