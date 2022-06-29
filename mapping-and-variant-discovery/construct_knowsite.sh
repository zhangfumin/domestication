#Conducting multi-sample variant calling using the 130 samples with high sequencing-depth
java -jar -Xmx100g -Xms100g -XX:+UseSerialGC GenomeAnalysisTK.jar -glm  SNP  -R rgp5.fa  -T  UnifiedGenotyper  -I <sample1.bam> -I <sample2.bam> .. -I <sample130.bam> -stand_call_conf 30 -stand_emit_conf 20   -o 130_samples.vcf

#Filtering the SNPs 
java -jar  GenomeAnalysisTK.jar  -R rgp5.fa  -T VariantFiltration  -V  130_samples.vcf  --filterExpression "FS >10, QD < 10, MQ< 40, AN<156, ReadPosRankSum <-0.5, BaseQRankSum<-0.5, MQRankSum<-0.5" --filterName "filter" -o knowsite.vcf
