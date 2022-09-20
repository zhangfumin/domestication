#Each individual ID (as defined in the VCF headerline) should be included on a separate line in the population_list file.
PopLDdecay -InVCF  <sample>.vcf -SubPop popolation_list  -MaxDist 300    -OutStat output_file
