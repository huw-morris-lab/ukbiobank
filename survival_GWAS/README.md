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

Run UKB_QC_script.txt in plink/R. Upload to Michigan Imputation server or similar.

Run postimput_bcftools_script.sh (I used kronos HPC).

Run postimput_plink_script.sh (I used kronos HPC).


# 3. Prepare for GWAS

Run preGWAS_script.txt in plink/R. This creates subsets of the genetic dataframe so that the GWAS can run in parallel.


# 4. Run GWAS

In parallel using kronos HPC.

Run gwas_make_rscripts.sh. This makes the number of R scripts according to the number of subsets of the data.

Run gwas_make_master_qsub.sh. This makes a single sh file with each lie as the qsub command which runs each R script.

Run master_qsub.sh. This submits all the jobs to kronos HPC.
