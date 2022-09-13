#calculating SNP allele frequence of populations
perl cal_pop_fr.pl  sam_info samples.vcf groups.fr

#groups.fr can divided into 12 frequence fils of 12 chromosomes   
#Building index for frequence files
perl build.fr.idx.pl group_chr01.fr

#calculating population parameters, pi, dxy and fst  
perl cal_pop_paramater.pl group_chr01.fre genes.info group_chr01_dis

#Constructing neighbor-join tree for each gene 
perl cal_pop_tree.pl group_chr01_dis group_chr01_dxy_tree group_chr01_fst_tree
 
#Performing 10,000 times of bootstrap resampling and generating 10,000 NJ trees
perl cal_bootstrip.pl gene_info 1000 1000 10000

