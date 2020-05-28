# ukbiobank
Pipeline for processing UK biobank clinical data

Date: May 2020

Last updated: 28/05/2020

Authors: Manuela Tan

## General description and purpose

This covers cleaning UK Biobank clinical data and pulling PD cases. 


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
