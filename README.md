# ukbiobank
Pipeline for processing UK biobank clinical data

Date: May 2020

Last updated: 01/11/2020

Authors: Manuela Tan

## General description and purpose

This covers cleaning UK Biobank clinical data and pulling PD cases. This is NOT a guide on how to download the data, detailed instructions have already been provided by UK Biobank. Note that these links and utilities may change/update over time. Refer to the email which says your data is available.

Data showcase to find out about data items that are available: https://biobank.ndph.ox.ac.uk/ukb/index.cgi
How to download UKB data: https://biobank.ctsu.ox.ac.uk/crystal/exinfo.cgi?src=AccessingData
How to download genetic data (also linked in the guide above): 
using ukbgene (older) - https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=664
using gfetch (newer version of download client) - https://biobank.ctsu.ox.ac.uk/showcase/refer.cgi?id=668

Once you have downloaded all the data, follow the steps below (as a guide).


# 1. Convert dataset

Convert dataset from .enc_ukb format to a workable format for R. Replace the file name with whatever your key is.

```
./ukbconv ukb23456.enc_ukb csv
```
or

```
./ukbconv ukb23456.enc_ukb r
```
The r format makes a .tab file and then a separate R script to recode all the categorical variables.

If making a csv file, note that the Data-Codings will be retained rather than replaced by their meanings.


Make data dictionary.
```
./ukbconv ukb23456.enc_ukb docs
```


# 2. Select relevant fields from the clinical dataset

As the dataset is very large and has lots of variables which I do not use, I run this R script to just save relevant variables. This outputs an R data format file .rds.

```
#!/usr/bin/Rscript

library(data.table)
library(dplyr)

#---Load in data---####

#Read in csv file
data <- fread("ukb37485.csv", header = T)


#---Select relevant columns---####

export <- data %>%
	select(eid, 
		`31-0.0`, #sex
		`34-0.0`, #year of birth
		`52-0.0`, #month of birth
		starts_with("53"), #date of baseline assessment
		starts_with("40023"), #DEATH and DEATH_CAUSE table. count of the number of rows for each participant in the DEATH table.
		#The DEATH table replicates the information in Field 40000, Field 40018 and Field 40020. The DEATH_CAUSE table replicates the information in Field 40001 and Field 40002
		#SUPPOSEDLY according to UK Biobank but the Data Portal tables are updated more frequently than the main dataset 
		starts_with("40000"), #date of death
		starts_with("40007"), #age at death
		starts_with("40001"), #cause of death primary
		starts_with("40014"), #cause of death primary - addendum
		starts_with("40002"), #secondary/contributing cause of death
		starts_with("40015"), #secondary/contributing cause of death addendum
		starts_with("40010"), #description of cause of death
		starts_with("40018"), #death record format
		starts_with("40020"), #death record origin
		starts_with("41270"), #HES summary diagnosis in ICD10. Summary of the distinct diagnosis codes both primary and secondary
		starts_with("41271"), #HES summary diagnosis in ICD9. Summary of the distinct diagnosis codes both primary and secondary
		starts_with("41202"), #Diagnoses - primary ICD10 (should be incorporated into 41270)
		starts_with("41204"), #Diagnoses - secondary ICD10 (should be incorporated into 41270)
		starts_with("41203"), #Diagnosis - primary in ICD9 (should be incorporated into 41271)
		starts_with("41205"), #Diagnosis - secondary in ICD9 (should be incorporated into 41271)
		starts_with("41280"), #date of first in-patient diagnosis
		starts_with("41262"), #Date of first in-patient diagnosis - main ICD10
		starts_with("41281"), #Date of first in-patient diagnosis - ICD9
		starts_with("20002"), #self-report non-cancer illness
		starts_with("20009"), #participant age for sefl-reported non-cancer illness
		starts_with("20008"), #year of the self-reported non-cancer illness
		starts_with("20013"), #date or age was recorded for self-reported illness
		starts_with("42032"), #Date of parkinson's disease report
		starts_with("42033"), #Source of PD report
		starts_with("42031"), #Source of all cause parkinsonism report
		starts_with("131023"), #Source of report of G20 (PD)
		starts_with("131022"), #Date G20 first reported (PD)
		)

#Save as R data file for analysis
saveRDS(export, "ukb_selectedFields.rds")
```

This takes a while (depending on how many variables you are interested in) so I run it using the HPC on kronos. 
```
qsub -pe make 2 -cwd ukb_selectFields_script_2020_06_13.R
```

# 2. Download death data from the Data Portal
According to UK Biobank, this data should be replicated in the main dataset but the Data Portal tables are updated more regularly than the main clinical dataset (Showcase Update only two or three times a year).

Further information on how to download mortality data is available in the UK Biobank access guides
https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/DeathLinkage.pdf

To access the Data Portal
* Access your project application through the AMS
* In the Downloads tab (where you downloaded all the other utilities to get the main dataset), go to the Data Portal tab 
* Connect to Record Repository

You can use SQL to query the tables, or just download the entire tables.

To get the full unique links, go to the Table Download tab and enter the name of the table you want. This will generate the wget commands like the ones below.

Download DEATH table
```
wget -nd -Odeath_cause.txt https://biota.ndph.ox.ac.uk/XXXXXX
```

Download DEATH_CAUSE table
```
wget -nd -Odeath.txt https://biota.ndph.ox.ac.uk/XXXXXX
```

