#Mapping short reads to reference Nipponbare genome (IRGSP build 5.0)
#Reference genome sequence file: rgp5.fa    

#Take whole genome resequencing data of sample nMYA_80724_2 as an example to demonstrate the mapping workflow
#Fastq file of raw reads for nMYA_80724_2: FCHYKJWCCXY_L7_wHAXPI104954-93_1.fq.gz FCHYKJWCCXY_L7_wHAXPI104954-93_2.fq.gz
#                                          FCHYKJWCCXY_L6_wHAXPI104954-93_1.fq.gz FCHYKJWCCXY_L6_wHAXPI104954-93_2.fq.gz

#Acquiring clean pire-end reads form raw data
trimmomatic PE -phred33 FCHYKJWCCXY_L7_wHAXPI104954-93_1.fq.gz FCHYKJWCCXY_L7_wHAXPI104954-93_2.fq.gz Clean_FCHYKJWCCXY_L7_wHAXPI104954-93_1.fq.gz Clean_FCHYKJWCCXY_L7_wHAXPI104954-93_1_unpaired.fq.gz Clean_FCHYKJWCCXY_L7_wHAXPI104954-93_2.fq.gz Clean_FCHYKJWCCXY_L7_wHAXPI104954-93_2_unpaired.fq.gz
trimmomatic PE -phred33 FCHYKJWCCXY_L6_wHAXPI104954-93_1.fq.gz FCHYKJWCCXY_L6_wHAXPI104954-93_2.fq.gz Clean_FCHYKJWCCXY_L6_wHAXPI104954-93_1.fq.gz Clean_FCHYKJWCCXY_L6_wHAXPI104954-93_1_unpaired.fq.gz Clean_FCHYKJWCCXY_L6_wHAXPI104954-93_2.fq.gz Clean_FCHYKJWCCXY_L6_wHAXPI104954-93_2_unpaired.fq.gz

#Aligning the short reads of the sample to the reference genome
bwa   mem  -aM    -R      '@RG\tID:L7\tSM:nMYA_80724_2\tPL:ILLUMINA\tLB:nMYA_80724_2\tPU:wHAXPI104954-93'   rgp5.fa     Clean_FCHYKJWCCXY_L7_wHAXPI104954-93_1.fq.gz  Clean_FCHYKJWCCXY_L7_wHAXPI104954-93_2.fq.gz > nMYA_80724_2_1.sam
bwa   mem  -aM    -R      '@RG\tID:L6\tSM:nMYA_80724_2\tPL:ILLUMINA\tLB:nMYA_80724_2\tPU:wHAXPI104954-93'   rgp5.fa     Clean_FCHYKJWCCXY_L6_wHAXPI104954-93_1.fq.gz  Clean_FCHYKJWCCXY_L6_wHAXPI104954-93_2.fq.gz > nMYA_80724_2_2.sam

#Merging sam files into one bam files
java -jar  /picard-tools-1.119/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT    INPUT=nMYA_80724_2_1.sam    INPUT=nMYA_80724_2_2.sam    OUTPUT=nMYA_80724_2.bam

#Marking duplicate reads
java -jar  /picard-tools-1.119/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES= false MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT=nMYA_80724_2.bam  OUTPUT=nMYA_80724_2_dedup.bam  METRICS_FILE=nMYA_80724_2_dedup.metrics

#Building index 
java -jar  /picard-tools-1.1119/BuildBamIndex.jar VALIDATION_STRINGENCY=LENIENT  I=nMYA_80724_2_dedup.bam

#Local realignment around indels
java -jar  GenomeAnalysisTK.jar  -R rgp5.fa   -T RealignerTargetCreator    -I  nMYA_80724_2_dedup.bam  -o nMYA_80724_2_dedup.intervals
java -jar  GenomeAnalysisTK.jar  -R rgp5.fa   -T IndelRealigner     -filterNoBases  -targetIntervals  nMYA_80724_2_dedup.intervals  -I nMYA_80724_2_dedup.bam  -o  nMYA_80724_2_dedup_realn.bam

#Base Quality Score Recalibration
#The knownsites file was obtained according to the script construct_knowsites.sh 
java -jar  GenomeAnalysisTK.jar  -R rgp5.fa  -T BaseRecalibrator -knownSites 130_wild_samples_knowsites.vcf   -I  nMYA_80724_2_dedup_realn.bam  -o  nMYA_80724_2_BQSR.grp
java -jar  GenomeAnalysisTK.jar  -R rgp5.fa  -T PrintReads    -I  nMYA_80724_2_dedup_realn.bam  -BQSR nMYA_80724_2_BQSR.grp -o nMYA_80724_2_BQSR.bam

