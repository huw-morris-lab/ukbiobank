#/bin/bash
plink --vcf UKB_PD.Michigan_HRC.allchromosomes.converted.R2_0.8.vcf.gz \
--double-id \
--allow-extra-chr 0 \
--maf 0.01 \
--make-bed \
--out UKB_PD.Michigan_HRC.allchromosomes.converted.R2_0.8.MAF_0.01