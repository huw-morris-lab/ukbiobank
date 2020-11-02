# PD mortality GWAS
Scripts for running PD progression GWAS of mortality in UKB.

Date: Nov 2020

Last updated: 01/11/2020

Authors: Manuela Tan

# 1. Process clinical data

Run R script ukb_clinical.R

This creates event and time to event variables and exports files with clinical data for analysis.

# 2. QC and imputation of genetic data

I performed my own QC and imputaton of UKB genotype data, rather than using the existing imputed data. This is because the PD participants may be from a slightly different population from the other UKB participants.

# 3. Prepare for GWAS

Run preGWAS_script.txt in plink/R.

# 4. Run GWAS

In parallel using kronos HPC
