#Discovering variants by applying a minimum base quality score of 20 and generating an intermediate GVCF file for each sample
java -jar  GenomeAnalysisTK.jar -R rgp5.fa  -T HaplotypeCaller -I <sample.bam> --emitRefConfidence GVCF  --min-base-quality-score 20 -maxAltAlleles 12  -variant_index_type LINEAR -variant_index_parameter 128000 -o <sample.vcf>

#Combining GVCF files of all samples into different groups of GVCF files
java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar   -R rgp5.fa  -T CombineGVCFs  --variant <sample_1.g.vcf> --variant <sample_2.g.vcf>  .. <sample_n.g.vcf>  --variant <sample_n.g.vcf> -o <combine_1.g.vcf>

#Discovering variants and genotype data of each chromosome
for i in {1..12}; do java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar  -U LENIENT_VCF_PROCESSING   --disable_auto_index_creation_and_locking_when_reading_rods  -R rgp5.fa -nt 3  -T GenotypeGVCFs  -L chr${i}build05r1.fasta  --variant <combine_1.g.vcf> --variant <combine_2.g.vcf> .. --variant <combine_n.g.vcf>  -o 1578_samples_hp_chr${i}.vcf; done

#Combining  vcf files of 12 chromosomes
java -jar -Xmx30g -Xms30g -XX:+UseSerialGC  picard-tools-1.119/MergeVcfs.jar  I=1578_samples_hp_chr01.vcf I=1578_samples_hp_chr02.vcf I=1578_samples_hp_chr03.vcf I=1578_samples_hp_chr04.vcf I=1578_samples_hp_chr05.vcf I=1578_samples_hp_chr06.vcf I=1578_samples_hp_chr07.vcf I=1578_samples_hp_chr08.vcf I=1578_samples_hp_chr09.vcf I=1578_samples_hp_chr10.vcf I=1578_samples_hp_chr11.vcf I=1578_samples_hp_chr12.vcf O=1578_samples_hp.vcf

#Selecting SNP
java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar  -R rgp5.fa -nt 8  -T SelectVariants -V 1578_samples_hp.vcf -selectType SNP -o 1578_samples_hp_SNP.vcf

#Variant Quality Score Recalibration <VQSR>
java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar   -R rgp5.fa  -T VariantRecalibrator  -nt 8  -input  1578_samples_hp_SNP.vcf  -resource:snp,known=true,training=true,truth=true,prior=10.0 knowsites.vcf -an DP  -an QD -an MQ  -an FS -an BaseQRankSum -an MQRankSum  -an ReadPosRankSum -mode SNP -tranche 99.9 -tranche 99.0     -tranche 97.0 -tranche 95.0 -tranche 93.0 -tranche 90.0  -recalFile 1578_samples_hp_SNP_VQSR.recal  -tranchesFile 1578_samples_hp_SNP_VQSR.tranches  -rscriptFile 1578_samples_hp_SNP_VQSR.plot.R
java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar   -R rgp5.fa  -T  ApplyRecalibration  -nt 8  -input  1578_samples_hp_SNP.vcf  -mode SNP  --ts_filter_level 90.0   -recalFile  1578_samples_hp_SNP_VQSR.recal  -tranchesFile 1578_samples_hp_SNP_VQSR.tranches  -o 1578_samples_hp_SNP_VQSR_90.vcf
#Filtering the SNPs with AN < 1578 (i.e., with less than half of samples with genotypes) 
java -jar  GenomeAnalysisTK.jar   -R rgp5.fa  -T VariantFiltration  -V  1578_samples_hp_SNP_VQSR_90.vcf  --filterExpression "AN < 1578" --filterName "filter" -o 1578_samples_snp.vcf

#Selecting INDEL and filtering INDEL
java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar   -R rgp5.fa -nt 8  -T SelectVariants -V 1578_samples_hp.vcf -selectType INDEL -o 1578_samples_hp_INDEL.vcf
java -jar  GenomeAnalysisTK.jar   -R rgp5.fa  -T VariantFiltration  -V  1578_samples_hp_INDEL.vcf  --filterExpression "QD < 5.0 || FS > 30.0 || MQ < 20.0 " --filterName "filter" -o 1578_samples_INDEL.vcf

