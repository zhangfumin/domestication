#date_200827_ind1_vcf2smc   
#wangmeixia

export PATH=/software/miniconda3/bin:$PATH
source activate /software/miniconda3/
conda activate smcpp

#take the ind1 population for example
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr01.gz ind_623.vcf.gz ind_chr01_input chr01build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357  # * represents the omitted individual number  
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr02.gz ind_623.vcf.gz ind_chr02_input chr02build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr03.gz ind_623.vcf.gz ind_chr03_input chr03build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr04.gz ind_623.vcf.gz ind_chr04_input chr04build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr05.gz ind_623.vcf.gz ind_chr05_input chr05build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr06.gz ind_623.vcf.gz ind_chr06_input chr06build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr07.gz ind_623.vcf.gz ind_chr07_input chr07build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr08.gz ind_623.vcf.gz ind_chr08_input chr08build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr09.gz ind_623.vcf.gz ind_chr09_input chr09build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr10.gz ind_623.vcf.gz ind_chr10_input chr10build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr11.gz ind_623.vcf.gz ind_chr11_input chr11build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
smc++ vcf2smc -m ind_623_SNPdensity_rm_chr12.gz ind_623.vcf.gz ind_chr12_input chr12build05r1.fasta ind1:ind_MAS_13910,ind_MAS_71644,ind_LAO_13010,******,ind_CAM_81357
                                                                                                                                                            
