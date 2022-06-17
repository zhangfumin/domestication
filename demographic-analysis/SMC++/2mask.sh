#date_ 200819_mask
#wangmeixia

#The loci with 0 SNPs in a window are masked
for ID in  "aro_28" "aus_97" "ind_623" "ray_13" "tej_101" "trj_181" "NIV1_93" "NIV2_86" "RUF1_87" "RUF2_142"        
do
echo $ID 
/build/bin/vcftools --vcf  "$ID".vcf --SNPdensity 10000 --out "$ID"_SNPdensity 
cat "$ID"_SNPdensity.snpden|awk '{if($3==0)print$1"\t"$2+1"\t"$2+10000}' > "$ID"_SNPdensity_cut.snpden
grep 'chr01build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr01
grep 'chr02build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr02
grep 'chr03build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr03
grep 'chr04build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr04
grep 'chr05build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr05
grep 'chr06build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr06
grep 'chr07build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr07
grep 'chr08build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr08
grep 'chr09build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr09
grep 'chr10build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr10
grep 'chr11build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr11
grep 'chr12build05r1.fasta' "$ID"_SNPdensity_cut.snpden > "$ID"_SNPdensity_rm_chr12
done

#compress and index the vcf files
for CHR in "chr01" "chr02" "chr03" "chr04" "chr05" "chr06" "chr07" "chr08" "chr09" "chr10" "chr11" "chr12"         
do
echo $CHR 
/build/bin/bgzip aro_28_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed aro_28_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip aus_97_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed aus_97_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip ind_623_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed ind_623_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip ray_13_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed ray_13_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip tej_101_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed tej_101_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip trj_181_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed trj_181_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip NIV1_93_SNPdensity_rm_"$CHR"            
/build/bin/tabix -p bed NIV1_93_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip NIV2_86_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed NIV2_86_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip RUF1_87_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed RUF1_87_SNPdensity_rm_"$CHR".gz
/build/bin/bgzip RUF2_142_SNPdensity_rm_"$CHR"             
/build/bin/tabix -p bed RUF2_142_SNPdensity_rm_"$CHR".gz
done
