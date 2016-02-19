Commands used for the VCF manipulation in Genome_exome files

./vcftools --vcf indels.recaliberated.vcf --remove-filtered-all --recode --recode-INFO-all --out Pass

## give 131404 out of 179005 sites

./vcftools --vcf Pass.recode.vcf --remove-indels --recode --recode_INFO-all --out SNP_only_Gen_Exo

## give 123847 SNPs out of 131404 sites

./vcftools --vcf Pass.recode.vcf --keep-only-indels --recode --recode_INFO-all --out SNP_only_Gen_Exo

## 10257 out of possible 131404 sites






Exome files only

./vcftools --vcf indels.recaliberated.vcf --remove-filtered-all --recode --recode-INFO-all --out Pass

## give us 134445 out of 157245

./vcftools --vcf Pass.recode.vcf --remove-indels --recode --recode_INFO-all --out SNP_only_Exo

## give us 123990 SNPs out of 134445 sites

./vcftools --vcf Pass.recode.vcf --keep-only-indels --recode --recode_INFO-all --out Indel_only_Exo

## give us 10455 out of possible 134445 sites

