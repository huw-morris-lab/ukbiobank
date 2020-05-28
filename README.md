# ukbiobank
Pipeline for processing UK biobank clinical data


1. Convert dataset

Convert dataset from .enc_ukb format to a workable format for R. Replace the file name with whatever your key is.

```
./ukbconv ukb23456.enc_ukb csv
```
or

```
./ukbconv ukb23456.enc_ukb r
```
The r format makes a .tab file and then a separate R script to code 


Make data dictionary.
```
./ukbconv ukb23456.enc_ukb docs
```
