#Output the genotype data in PLINK PED format
vcftools --vcf  1493_samples.vcf   --plink-tped   --out 1493_samples
#Creating a binary ped file
plink --tfile 1493_samples --recode --out 1493_samples
plink  --file 1493_samples    --make-bed --out 1493_samples_pruned --recode
#Principal component analysis. Get the GRM(1493_samples_pruned_grm) and output the first 20 (--pca 20) eigenvectors and all the eigenvalues.
gcta64 --bfile  1493_samples_pruned  --make-grm  --out 1493_samples_pruned_grm
gcta64  --grm 1493_samples_pruned_grm  --pca 20 --out 1493_samples_pca_20
#Estimating the pairwise genetic distance matrix
plink --bfile  1493_samples_pruned  --distance 1-ibs flat-missing --out  1493_samples_pruned_distance
#removing each SNP that has an R2 value of greater than 0.1 with any other SNP within a 500-SNP sliding window(step 200-SNP)
plink --noweb  --file 1493_samples    --indep-pairwise 500  200 0.1
#copying the remaining (untargetted) SNPs to 1493_samples_thin_pruned.bed file
plink --noweb --file 1493_samples   --extract plink.prune.in --make-bed --out 1493_samples_thin_pruned
#running ADMIXTURE with cross-validation for K values 1 to 12
for K in {1..12};admixture --cv 1493_samples_thin_pruned.bed ${K} > 1493_samples_log${K}.out














