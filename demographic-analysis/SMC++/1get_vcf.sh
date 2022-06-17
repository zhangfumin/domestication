#date_20200817_get_vcf
#wangmeixia

#get_vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep aus_97  > aus_97.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep aro_28  > aro_28.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep ray_13  > ray_13.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep tej_101 > tej_101.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep trj_181 > trj_181.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep ind_623 > ind_623.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep NIV1_93 > NIV1_93.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep NIV2_86 > NIV2_86.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep RUF1_87 > RUF1_87.vcf
/build/bin/vcftools --vcf  200731_1586_sample_hp_SNP_VQSR_90_filter_pass.vcf --recode --recode-INFO-all --stdout --keep RUF2_142> RUF2_142.vcf

#compress and index the vcf files
for ID in  "aro_28" "aus_97" "ind_623" "ray_13" "tej_101" "trj_181" "NIV1_93" "NIV2_86" "RUF1_87" "RUF2_142"        
do
echo $ID 
/build/bin/bgzip -c "$ID".vcf >  "$ID".vcf.gz    
/build/bin/tabix -p vcf "$ID".vcf.gz             
done
