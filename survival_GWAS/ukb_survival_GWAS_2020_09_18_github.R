#------------SURVIVAL GWAS IN UK BIOBANK------------#
#Using this documentation to determine PD cases and incident vs. prevalent
#https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/alg_outcome_pdp.pdf

#---Load packages---####
library(ggplot2)
library(readstata13)
library(gridExtra)
library(gtable)
library(grid)
library(reshape2)
library(pROC)
library(ROCR)
library(readxl)
library(survival)
library(survminer)
library(lubridate)
library(lme4)
library(car)
library(MASS)
library(factoextra)
library(data.table)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(corrplot)
library(VennDiagram)

#---Load functions---####

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}


## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}

#Normalise function
normFunc <- function(x){
  (x-mean(x, na.rm = T))/sd(x, na.rm = T)
}

#---Load clinical data---####

#Selected fields from the full dataset
data <- readRDS(file = "../ukb_selectedFields.rds")

#FIELD NAMES
# `31-0.0`, #sex
# `34-0.0`, #year of birth
# `52-0.0`, #month of birth
# starts_with("53"), #date of baseline assessment
# starts_with("40023"), #DEATH and DEATH_CAUSE table. count of the number of rows for each participant in the DEATH table.
# The DEATH table replicates the information in Field 40000, Field 40018 and Field 40020. The DEATH_CAUSE table replicates the information in Field 40001 and Field 40002
# starts_with("40000"), #date of death
# starts_with("40007"), #age at death
# starts_with("40001"), #cause of death primary
# starts_with("40014"), #cause of death primary - addendum (NOT IN DATA)
# starts_with("40002"), #secondary/contributing cause of death
# starts_with("40015"), #secondary/contributing cause of death addendum (NOT IN DATA)
# starts_with("40010"), #description of cause of death
# starts_with("40018"), #death record format
# starts_with("40020"), #death record origin
# starts_with("41270"), #HES summary diagnosis in ICD10. Summary of the distinct diagnosis codes both primary and secondary
# starts_with("41271"), #HES summary diagnosis in ICD9. Summary of the distinct diagnosis codes both primary and secondary
# starts_with("41202"), #Diagnoses - primary ICD10 (should be incorporated into 41270)
# starts_with("41204"), #Diagnoses - secondary ICD10 (should be incorporated into 41270)
# starts_with("41203"), #Diagnosis - primary in ICD9 (should be incorporated into 41271)
# starts_with("41205"), #Diagnosis - secondary in ICD9 (should be incorporated into 41271)
# starts_with("41280"), #date of first in-patient diagnosis. Date of 41270
# starts_with("41262"), #Date of first in-patient diagnosis - main ICD10. Date of 41202
# starts_with("41281"), #Date of first in-patient diagnosis - ICD9. Date of first diagnosis of 41203
# starts_with("20002"), #self-report non-cancer illness
# starts_with("20009"), #participant age for sefl-reported non-cancer illness
# starts_with("20008"), #year of the self-reported non-cancer illness
# starts_with("20013"), #date or age was recorded for self-reported illness
# starts_with("42032"), #Date of parkinson's disease report
# starts_with("42033"), #Source of PD report
# starts_with("42031") #Source of all cause parkinsonism report
# starts_with("131023"), #Source of report of G20 (PD)
# starts_with("131022"), #Date G20 first reported (PD)

#---Read in ICD codes---####

ICD_10 <- fread("../coding19.tsv", header = T)

ICD_9 <- fread("../coding87.tsv", header = T)

#---Read in death data separate from main dataset---#####
#This is from the Data Portal: Record Repository
#Not sure if this death data is also in the main database - need to merge and check

death <- fread("../death.txt")
death_cause <- fread("../death_cause.txt")

#---Save and load RDS---####

#Make filename with today's date
filename_rds <- paste("workspace_",today(), ".Rdata", sep = "")

#Save whole workspace
save.image(filename_rds)

#Load workspace
load("workspace_2020-06-17.Rdata")

#---Select PD cases from HES---####
#Defined using HES, self-report, or cause of death
#I is the instance index and A the array index

#Cases from hospital episode statistics
#This is a summary across all timepoints from ICD10 codes
PD_HES_ICD10 <- data %>% 
  filter_at(vars(contains("41270")), any_vars(.=="G20"))
#Note that this variable 41270 should incorporate data from 41202 (main ICD10 diagnosis) and 41204 (secondary ICD10 diagnosis)

#Check that the numbers match using the original ICD10 primary and secondary variables
PD_HES_ICD10_check <- data %>% 
  filter_at(vars(contains("41202"), contains("41204")), any_vars(.=="G20"))

#Cases from hospital episode statistics from ICD9 codes
#This is a summary across all timepoints from ICD9 codes
PD_HES_ICD9 <- data %>% 
  filter_at(vars(contains("41271")), any_vars(.=="332"))

#---Select PD cases from self-report---####

#Cases from self-report
#PD code is 1262
#This includes people who also have the ICD diagnosis
PD_selfreport <- data %>%
  filter_at(vars(contains("20002")), any_vars(.=="1262"))

#---Select PD cases from death data---####

#Cases from death data - using both primary and secondary causes
PD_death <- data %>% 
  filter_at(vars(contains("40001"), contains("40002")), any_vars(.=="G20"))

#Check freetext death field
# check <- PD_death %>% 
#   select(contains("40010"), contains("40001"), contains("40002"))

#Check to see if there are any cases with PD coded in the freetext field that are not in the ICD codes
freetext_death <- data %>% 
  filter(str_detect(`40010-0.0`, "Parkinson")) %>% 
  select(eid, "40010-0.0", contains("40001"), contains("40002"))

PDdeath_missing <- freetext_death %>% 
  anti_join(PD_death, by = "eid")
#This looks fine - not PD cases

#---Check against death and death cause tables---####
#Comparing death data in main table to the separate death and death cause tables

#Filter death cause table for G20
death_cause_G20 <- death_cause %>% 
  filter(cause_icd10 == "G20")

#Check there are no duplicated IDs
death_cause_G20$eid[duplicated(death_cause_G20$eid)]

#The death cause tables have more individuals than the PD_death (40001 and 40002) in the main data
#These are new tables - possibly more updated

#Get date of death for the PD patients who have died
death_cause_G20_date <- death_cause_G20 %>% 
  left_join(death, by = "eid")

#There is one duplicated record 4161184 but the date of death is the same
death_cause_G20_date$eid[duplicated(death_cause_G20_date$eid)]

#Keep just unique IDs
death_cause_G20_date_unique <- death_cause_G20_date %>% 
  distinct(eid, .keep_all = TRUE)

#Why are there cases in the death cause tables that are not in the main data?
death_cause_G20_date_unique %>% 
  anti_join(PD_death, by = "eid")
#Some of these could be HES PD cases who have died

#---Convert death tables into wide format---####
#This will make it easier to merge with the main dataset
#So there is only one row per participant

#First check for duplicated records in death table
death$eid[duplicated(death$eid)]
#There are a handful of duplicated records

#Count how many individuals have a second death record
death %>% 
  group_by(ins_index) %>% 
  summarise(count = n())

#Arrange by ID and ins_index and then keep only unique records (only ins_index==0)
death_unique <- death %>% 
  arrange(eid, ins_index) %>% 
  distinct(eid, .keep_all = TRUE)

#Convert death_cause table from long to wide
death_cause_wide <- death_cause %>% 
  select(-level) %>% 
  spread(arr_index, cause_icd10)
#The primary cause of death (level 1) is always in the arr_index 0 position
#So don't really need to keep the level variable

death_cause_wide %>% 
  group_by(ins_index) %>% 
  summarise(count = n())

#Rename all the cause of death variables (0 to 14)
#Also keep only the records with ins_index == 0 
#The ins_index is to indicate duplicated death records and don't know which is right when they differ
death_cause_wide_unique <- death_cause_wide %>% 
  rename_if(is.character, .funs = funs(paste0("cause_", .))) %>% 
  filter(ins_index == 0)

#Merge death dates and death cause data
death_dates_causes <- death_unique %>%
  left_join(death_cause_wide_unique, by = c("eid", "ins_index"))

#There are a handful of individuals who are missing primary and any other cause of death
#But still want to include these as we have the date of death
death_dates_causes %>% 
  filter(is.na(cause_0))

#---Calculate date of PD diagnosis from PD HES ICD10 diagnosis---####
#Date of first HES diagnosis ICD10 - date of first diagnosis of 41270

#Select just the HES data - diagnosis and date
PD_HES_ICD10_selected <- PD_HES_ICD10 %>% 
  select(eid, contains("41270"), contains("41280"))

#Get the variable names for the ICD10 diagnosis and dates
ICD10_names <- colnames(PD_HES_ICD10_selected %>% 
                          select(contains("41270")))

ICD10_date_names <- colnames(PD_HES_ICD10_selected %>% 
                               select(contains("41280")))

#Convert to dataframe
PD_HES_ICD10_selected_df <- as.data.frame(PD_HES_ICD10_selected)

#Reshape to long format
#Gather both ICD values and date values
PD_HES_ICD10_selected_df_long <- reshape(PD_HES_ICD10_selected_df, idvar="eid", direction="long", 
                                         varying = list(ICD10_names, ICD10_date_names),
                                         v.names = c("ICD10", "date"))

#Arrange by ID
PD_HES_ICD10_selected_df_long <- PD_HES_ICD10_selected_df_long %>% 
  arrange(eid)

#Keep just date corresponding to PD diagnosis
PD_HES_ICD10_selected_df_long_G20 <- PD_HES_ICD10_selected_df_long %>% 
  filter(ICD10 == "G20")

#Check that each person only has one row
PD_HES_ICD10_selected_df_long_G20$eid[duplicated(PD_HES_ICD10_selected_df_long_G20$eid)]

#Rename date variable
PD_HES_ICD10_selected_df_long_G20 <- PD_HES_ICD10_selected_df_long_G20 %>% 
  rename(ICD10_G20_date = date) %>% 
  select(eid, ICD10_G20_date)

#Select just eid and date to merge
PD_HES_ICD10_selected_df_long_G20_merge <- PD_HES_ICD10_selected_df_long_G20 %>% 
  select(eid, ICD10_G20_date)

#Merge with main PD HES dataset
PD_HES_ICD10 <- PD_HES_ICD10 %>% 
  left_join(PD_HES_ICD10_selected_df_long_G20_merge, by = "eid")

#---Calculate date of PD self-report---####
#Date of self-report PD
# starts_with("20009"), #participant age for sefl-reported non-cancer illness
# starts_with("20008"), #year of the self-reported non-cancer illness
# 20008 is the interpolated time when the participant indicated the corresponding condition was first diagnosed by a doctor, measured in years.
# starts_with("20013"), #date or age was recorded for self-reported illness


#Select just the self-report data - diagnosis and date fields
#There is both age and year available - select both
PD_selfreport_selected <- PD_selfreport %>% 
  select(eid, contains("20002"), contains("20009"), contains("20008"))

#Get the variable names for the self-report diagnosis and dates
selfdiagnosis_names <- colnames(PD_selfreport_selected %>% 
                                  select(contains("20002")))

selfdiagnosis_age_names <- colnames(PD_selfreport_selected %>% 
                                      select(contains("20009")))

selfdiagnosis_year_names <- colnames(PD_selfreport_selected %>% 
                                       select(contains("20008")))

#Convert to dataframe
PD_selfreport_selected_df <- as.data.frame(PD_selfreport_selected)

#Reshape to long format
#Gather both selfreport diagnosis values, year values, and age values
PD_selfreport_selected_df_long <- reshape(PD_selfreport_selected_df, idvar="eid", direction="long", 
                                          varying = list(selfdiagnosis_names, selfdiagnosis_age_names, selfdiagnosis_year_names),
                                          v.names = c("selfreport", "age", "year"))

#Code visit where illness was reported, arrange by ID
PD_selfreport_selected_df_long <- PD_selfreport_selected_df_long %>% 
  mutate(instance = ifelse(time >= 1 & time <=34, "BL",
                           ifelse(time >=35 & time <= 68, "Visit1",
                                  ifelse(time >=69 & time <= 102, "Visit2",
                                         ifelse(time >=103 & time <=136, "Visit3", NA))))) %>% 
  arrange(eid)

#Keep just date corresponding to PD diagnosis
PD_selfreport_selected_df_long_1262 <- PD_selfreport_selected_df_long %>% 
  filter(selfreport == "1262")

#Check that each person only has one row
PD_selfreport_selected_df_long_1262$eid[duplicated(PD_selfreport_selected_df_long_1262$eid)]

#There are some duplicates - just take the first instance which is the earliest self-report date
#Based on ALG document
PD_selfreport_selected_df_long_1262_unique <- PD_selfreport_selected_df_long_1262 %>% 
  distinct(eid, .keep_all = TRUE)

#Rename date variable
PD_selfreport_selected_df_long_1262_unique_merge <- PD_selfreport_selected_df_long_1262_unique %>% 
  rename(selfreport_1262_age = age,
         selfreport_1262_year = year,
         selfreport_1262_instance = instance) %>% 
  select(eid, selfreport, selfreport_1262_age, selfreport_1262_year, selfreport_1262_instance)

#Merge with main PD HES dataset
PD_selfreport <- PD_selfreport %>% 
  left_join(PD_selfreport_selected_df_long_1262_unique_merge, by = "eid")

#Calculate date from year decimal
PD_selfreport <- PD_selfreport %>% 
  mutate(selfreport_1262_date = format(date_decimal(selfreport_1262_year), "%Y-%m-%d")) #make date using the 1st of each mon


#---Subsets of PD cases defined by HES vs. self-report vs. death---####

#Identify PD cases from HES that are not in self report or death
PD_HES_ICD10_only <- PD_HES_ICD10 %>% 
  anti_join(death_cause_G20_date_unique, by = "eid") %>% 
  anti_join(PD_selfreport, by = "eid")

#Identify PD cases from death data that are not in HES or self report
PD_death_only <- death_cause_G20_date_unique %>% 
  anti_join(PD_HES_ICD10, by = "eid") %>% 
  anti_join(PD_selfreport, by = "eid")

#Identify PD cases from self-report data that are not in HES or death
PD_selfreport_only <- PD_selfreport %>% 
  anti_join(PD_HES_ICD10,  by = "eid") %>% 
  anti_join(death_cause_G20_date_unique, by = "eid")

#PD cases that have self-reported and are in HES but not death
PD_HES_selfreport <- PD_HES_ICD10 %>% 
  inner_join(PD_selfreport, by = "eid") %>% 
  anti_join(death_cause_G20_date_unique, by = "eid") %>% 
  select(eid)

#PD cases that have all 3 reports - HES, self-report and death
PD_HES_selfreport_death <- PD_HES_ICD10 %>% 
  inner_join(PD_selfreport, by = "eid") %>% 
  inner_join(death_cause_G20_date_unique, by = "eid") %>% 
  select(eid)

#PD cases that have HES and death only, not selfreport
PD_HES_death <- PD_HES_ICD10 %>% 
  inner_join(death_cause_G20_date_unique, by = "eid") %>% 
  anti_join(PD_selfreport, by = "eid") %>% 
  select(eid)

#PD cases that have both self-report and death only, not HES
PD_selfreport_death <- PD_selfreport %>% 
  inner_join(death_cause_G20_date_unique, by = "eid") %>% 
  anti_join(PD_HES_ICD10, by = "eid")


#---Venn diagram of where PD cases have been identified from HES vs. self-report vs. death---####

#Venn diagram plot
#Appearance taken from https://www.r-graph-gallery.com/14-venn-diagramm.html
venn.diagram(
  x = list(PD_HES_ICD10$eid, PD_selfreport$eid, death_cause_G20_date_unique$eid),
  category.names = c("HES" , "self-report " , "death"),
  filename = 'plots/UKB_PD.png',
  output = TRUE ,
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
  rotation = 1
)

#---Combine all datasets for all PD cases---####
#Make dataset with all PD cases

#For the cases that have both HES and self-report data, want to include both dates

#Select just the self-report dates from the PD self-report cases
PD_selfreport_dateOnly <- PD_selfreport %>% 
  select(eid, selfreport_1262_date, selfreport_1262_instance) %>% 
  mutate(selfreport_PD = "yes")

#Merge  PD HES dataset and self-report dates - this will cover people who have both HES and self-report
PD_HES_ICD10_selfreportDate <- PD_HES_ICD10 %>% 
  mutate(HES_PD = "yes") %>% 
  left_join(PD_selfreport_dateOnly, by = "eid")

#Select variables for merge from HES dataset
PD_HES_ICD10_selfreportDate_merge <- PD_HES_ICD10_selfreportDate %>% 
  left_join(death_dates_causes, by = "eid") %>% #merge with updated death tables data
  select(eid, `31-0.0`, #sex
         `34-0.0`, #year of birth
         `52-0.0`, #month of birth
         starts_with("53-"), #date of baseline assessment
         dsource, source, date_of_death,
         starts_with("cause_"), 
         HES_PD, ICD10_G20_date,
         selfreport_PD, selfreport_1262_date, selfreport_1262_instance)

#Select variables for merge from selfreport dataset
PD_selfreport_merge <- PD_selfreport %>% 
  #Remove the cases that also have HES PD because we have already included their data in the HES dataset
  anti_join(PD_HES_ICD10_selfreportDate_merge, by = "eid") %>% 
  mutate(selfreport_PD = "yes") %>% 
  left_join(death_dates_causes, by = "eid") %>% #merge with updated death tables data
  #Select the same variables as we selected from the HES dataset
  select(eid, `31-0.0`, #sex
         `34-0.0`, #year of birth
         `52-0.0`, #month of birth
         starts_with("53-"), #date of baseline assessment
         dsource, source, date_of_death,
         starts_with("cause_"),
         selfreport_PD, selfreport_1262_date, selfreport_1262_instance) %>% 
  #Make some blank variables for HES stats - these patients do not have
  mutate(HES_PD = NA, 
         ICD10_G20_date = NA)

#Select variables for merge from death dataset
#Using the list of individuals from the death tables, not the main database
#This is just cases that have PD reported in death but not HES and not selfreport
PD_death_only_merge <- PD_death_only %>%
  select(eid) %>% 
  left_join(data, by = "eid") %>%
  left_join(death_dates_causes, by = "eid") %>% 
  #Select the same variables as we selected from the HES dataset
  select(eid, `31-0.0`, #sex
         `34-0.0`, #year of birth
         `52-0.0`, #month of birth
         starts_with("53-"), #date of baseline assessment
         dsource, source, date_of_death,
         starts_with("cause_")) %>% 
  #Make some blank variables for HES stats - these patients do not have
  mutate(HES_PD = NA, 
         ICD10_G20_date = NA,
         selfreport_PD = NA, 
         selfreport_1262_date = NA,
         selfreport_1262_instance = NA)


#Combine all three datasets - there should not be any duplicates/overlapping cases
PD_merged <- rbind(PD_HES_ICD10_selfreportDate_merge, PD_selfreport_merge, PD_death_only_merge)

#Check there are no duplicated IDs
PD_merged$eid[duplicated(PD_merged$eid)]

#---Check how PD cases have been classified vs. algorithmically defined in UKB---####

#By my definition
PD_merged %>% 
  group_by(HES_PD, selfreport_PD) %>% 
  summarise(count = n())

#By UKB definition - field 42033 Source of parkinson's disease report
data %>% 
  group_by(`42033-0.0`) %>% 
  summarise(count = n())
#This looks close enough - I think this is an earlier version (last updated 05 Mar 2020)
#https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=42033
#So there are slightly more PD cases in the current dataset 
#Note that with the updated death data there is a bit more mismatch with the numbers - more PD cases identified with death

#---Check for missing dates of PD report (either HES or selfreport)---####

#Check if there are any missing dates for HES PD
PD_merged %>% 
  filter(is.na(ICD10_G20_date), HES_PD == "yes") %>% 
  select(eid, ICD10_G20_date, HES_PD)
#No missing dates for HES PD

#Check if there are any missing dates for self-report PD
PD_merged %>% 
  filter(is.na(selfreport_1262_date), selfreport_PD == "yes") %>% 
  select(eid, selfreport_1262_date, selfreport_PD)
#No missing dates for selfreport PD

#Check if there are any dates that are 1900 - this is the missing value
dates <- PD_merged %>% 
  select(ICD10_G20_date, selfreport_1262_date) %>% 
  filter(ICD10_G20_date < 1903 | selfreport_1262_date < 1903)

#---Check if there are PD cases defined in HES or self report that have died but PD is not a cause of death---####

PD_merged %>% 
  filter(!is.na(date_of_death)) %>% 
  summarise(count = n())

#Select PD cases that are in HES and/or self-report that have died and have G20 in the cause of death
PD_merged_HES_self_PDdeath <- PD_merged %>% 
  filter(HES_PD == "yes" | selfreport_PD == "yes") %>% 
  filter(!is.na(date_of_death)) %>% 
  filter_at(vars(contains("cause_")), any_vars(.=="G20"))
#440 cases that have PD in HES or self report that have died also have PD as a cause of death

test <- PD_merged %>% 
  select(contains("cause_"))

#Find PD cases that are in HES or self-report that have died and do NOT have G20 in the cause of death
PD_merged_HES_self_noPDdeath <- PD_merged %>% 
  filter(HES_PD == "yes" | selfreport_PD == "yes") %>% 
  filter(!is.na(date_of_death)) %>% 
  anti_join(PD_merged_HES_self_PDdeath, by = "eid")

#Find PD cases that are identified from death only (not reported in HES or self-report)
PD_merged_PDdeath_only <- PD_merged %>% 
  filter(is.na(HES_PD) & is.na(selfreport_PD))

#---Determine if PD cases are incident or prevalent---####
#All defined using https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/alg_outcome_pdp.pdf


#PREVALENT PD
#Prevalent PD is when the first ICD code date is prior to the date of baseline assessment
#OR The participant has self-reported the condition at the baseline assessment, but the first ICD code date is post the date of baseline assessment.
#OR Parkinson’s Disease by self-report only: The participant has self-reported Parkinson’s Disease at baseline assessment, but without evidence of Parkinson’s Disease from linked HES APC, SMR01 or PEDW data (as defined above).

#INCIDENT PD (everything else)
#Excluding those with Parkinson’s Disease detected prior to baseline assessment
#Parkinson’s Disease detected by hospital admission EHR ICD10 code date post the baseline assessment.
#OR Parkinson’s Disease detected by death register only: No ICD codes but one (or more) ICD codes in death register records, in the underlying cause or any other position.

#Note that if patients are identified by self report ONLY after the baseline visit, the ALG guidelines do not include these as either incident or prevalent
#May be because there is bias in follow-up
PD_merged <- PD_merged %>% 
  mutate(PD_prevalent = ifelse(selfreport_1262_instance!="BL" & is.na(ICD10_G20_date), "undefined",
                               ifelse(selfreport_1262_instance == "BL" & is.na(ICD10_G20_date), "prevalent", 
                                      ifelse(ICD10_G20_date < `53-0.0`, "prevalent",
                                             ifelse(selfreport_1262_instance == "BL" & ICD10_G20_date > `53-0.0`, "prevalent", NA))))) %>% 
  mutate(PD_status = ifelse(is.na(PD_prevalent), "incident", 
                            ifelse(PD_prevalent == "prevalent", "prevalent",
                                   ifelse(PD_prevalent == "undefined", "undefined", NA))))


#Summary of prevalent vs. incident vs. undefined cases
PD_merged %>% 
  group_by(PD_status) %>% 
  summarise(count = n())

#Write list of incident and prevalent cases to extract from genotyping data
PD_incident <- PD_merged %>%
  filter(PD_status == "incident") %>%
  mutate(IID = paste(eid, "_", eid, sep = ""),
         FID = IID) %>%
  select(FID, IID)

write.table(PD_incident, "outputs/incident_IDs.txt", quote = F, row.names = F, col.names = F)

PD_prevalent <- PD_merged %>%
  filter(PD_status == "prevalent") %>%
  mutate(IID = paste(eid, "_", eid, sep = ""),
         FID = IID) %>%
  select(FID, IID)

write.table(PD_prevalent, "outputs/prevalent_IDs.txt", quote = F, row.names = F, col.names = F)


#---Set dates for PD diagnosis---####
#According to ALG document

#For PREVALENT CASES
#If a participant has both an ICD code and a self-report code, the earliest recorded date regardless of source is used.
#If a participant has both an ICD code and a self-report code, but the self-reported date is missing, the ICD code date is used 
#If the participant has ICD code(s) only, the earliest ICD code date is used.
#If the participant has self-report code(s) only, the earliest self-reported date is used

#For INCIDENT CASES
#If a participant has ICD codes in both hospital admission and death register records, the earliest recorded code date regardless of source is used.
#If ICD code(s) recorded in hospital admission only, the earliest ICD code date is used.
#If ICD code(s) recorded in death register only, the date of death is used.

#Using the new date_of_death variable which is from the DEATH tables, not the main database (40000 variable)
PD_merged <- PD_merged %>%
  mutate(PD_date = ifelse(PD_status == "prevalent" & is.na(selfreport_1262_date) & !is.na(ICD10_G20_date), ICD10_G20_date,
                          ifelse(PD_status == "prevalent" & !is.na(selfreport_1262_date) & is.na(ICD10_G20_date), selfreport_1262_date,
                                 ifelse(PD_status == "prevalent" & ICD10_G20_date < selfreport_1262_date, ICD10_G20_date,
                                        ifelse(PD_status == "prevalent" & ICD10_G20_date > selfreport_1262_date, selfreport_1262_date,
                                               ifelse(PD_status == "incident" & is.na(date_of_death) & !is.na(ICD10_G20_date), ICD10_G20_date,
                                                      ifelse(PD_status == "incident" & !is.na(date_of_death) & is.na(ICD10_G20_date), date_of_death,
                                                             ifelse(PD_status == "incident" & !is.na(date_of_death) & !is.na(ICD10_G20_date), ICD10_G20_date, 
                                                                    ifelse(PD_status == "undefined", selfreport_1262_date, "check")))))))))



#Convert PD diagnosis date to a date
PD_merged <- PD_merged %>% 
  mutate(PD_date = as.Date(PD_date))

#---Define last follow-up date---####

#Convert date of death to a date
PD_merged <- PD_merged %>% 
  mutate(date_of_death = as.Date(date_of_death, format = "%d/%m/%Y"))

#Find the latest death date
test <- PD_merged %>% 
  filter(!is.na(date_of_death)) %>% 
  arrange(desc(date_of_death)) %>% 
  select(eid, date_of_death)
#The last date of death is 21/04/2020

#---Determine PD cause of death in incident vs. prevalent cases---####

#COunt how many incident PD cases died
PD_merged %>% 
  filter(PD_status == "incident") %>% 
  filter(!is.na(date_of_death)) %>% 
  summarise(count = n())

#Count how many incident PD cases died - excluding those 129 who were classified from death records only
PD_merged %>% 
  filter(PD_status == "incident") %>% 
  filter(!is.na(date_of_death)) %>% 
  filter(HES_PD == "yes" | selfreport_PD == "yes") %>%
  summarise(count = n())

#Count how many incident PD cases have PD as a cause of death
PD_merged %>% 
  filter(PD_status == "incident") %>% 
  filter(!is.na(date_of_death)) %>% 
  filter(HES_PD == "yes" | selfreport_PD == "yes") %>% 
  filter_at(vars(contains("cause_")), any_vars(.=="G20")) %>% 
  summarise(count = n())

#Count how many prevalent PD cases have PD as a cause of death
PD_merged %>% 
  filter(PD_status == "prevalent") %>% 
  filter(!is.na(date_of_death)) %>% 
  summarise(count = n())

#---Make new variable for PD-related death cases vs. interrupted---####

#Code people who have PD (G20) listed as one of the causes of death
PDdeath <- PD_merged %>% 
  filter(!is.na(date_of_death)) %>% 
  filter(cause_0 =="G20") %>% 
  mutate(PDdeath = "yes") %>% 
  select(eid, PDdeath)

PD_merged <- PD_merged %>% 
  left_join(PDdeath, by = "eid")
  
PD_merged %>% 
  filter(!is.na(date_of_death)) %>% 
  group_by(PD_status, PDdeath) %>% 
  summarise(count = n())

########## TABLE 1: SUMMARY DEMOGRAPHICS ########## 
#---Gender summaries---####

PD_merged <- PD_merged %>% 
  mutate(gender = ifelse(`31-0.0` == 0, "female",
                         ifelse(`31-0.0` == 1, "male", NA)))

PD_merged %>% 
  group_by(gender) %>% 
  summarise(count = n())

#Gender summaries grouped by PD status
PD_merged %>% 
  group_by(PD_status, gender) %>% 
  summarise(count = n())

#---Calculate age at PD diagnosis---####

#Calculate age at diagnosis
PD_merged <- PD_merged %>%
  mutate(dob = paste(`34-0.0`,`52-0.0`, "01", sep = "-"), #make DOB using the 1st of the month and year of birth
         dob = as.Date(dob)) %>% 
  mutate(age_diagnosis = dob %--% PD_date,
         age_diagnosis = as.duration(age_diagnosis) / dyears(1)) #calculate age at diagnosis in years

#Check individuals who are missing age at diagnosis
PD_merged %>% 
  filter(is.na(age_diagnosis)) %>% 
  summarise(count = n())
#These are the 129 PD cases identified from death records only

#Summary of PD age at diagnosis
PD_merged %>% 
  filter(!is.na(age_diagnosis)) %>% 
  summarise(mean_aad = mean(age_diagnosis, na.rm = TRUE),
            sd_aad = sd(age_diagnosis, na.rm = TRUE))

#Summary of PD age at diagnosis grouped by PD status
PD_merged %>% 
  filter(!is.na(age_diagnosis)) %>%
  group_by(PD_status) %>% 
  summarise(mean_aad = mean(age_diagnosis, na.rm = TRUE),
            sd_aad = sd(age_diagnosis, na.rm = TRUE))


#---Calculate disease duration at study entry for prevalent cases---####

#53-0.0 is the date of baseline assessment
PD_merged <- PD_merged %>% 
  mutate(date_entry = as.Date(`53-0.0`, format = "%Y-%m-%d")) %>% 
  mutate(age_entry = as.duration(dob %--% date_entry)/dyears(1)) %>% 
  mutate(disease_duration_onset = as.numeric(ifelse(PD_status == "incident", NA,
                                         ifelse(PD_status == "undefined", NA,
                                                ifelse(PD_status == "prevalent", 
                                                       (as.duration(PD_date %--% date_entry)/dyears(1)), "check")))))


########## MORTALITY AND SURVIVAL OUTCOMES ########## 
#---Create variables for event and time to event---####

#Create event and time to event variables
#The time to event is measured from PD diagnosis, not study entry
PD_merged <- PD_merged %>% 
  mutate(event_death = ifelse(!is.na(date_of_death), 1, 0)) %>% 
  #Get the interval between the PD date and date of death, or last follow-up date 01/04/2020
  mutate(timeToEvent_death = ifelse(!is.na(date_of_death), PD_date %--% date_of_death,
                                    ifelse(is.na(date_of_death), PD_date %--% as.Date("2020-04-01"), NA))) %>% 
  #Convert the interval to years
  mutate(timeToEvent_death = as.duration(timeToEvent_death) / dyears(1))

#Check individuals who are missing time to event data
check <- PD_merged %>% 
  select(eid, event_death, timeToEvent_death, date_of_death) %>% 
  filter(is.na(timeToEvent_death))
#These are the 129 people who are only identified as having PD in the death records and no other sources

#---Summary of mortality outcomes---####

#Count number of PD cases who have died overall (including those only identified from death records)
PD_merged %>% 
  group_by(event_death) %>% 
  summarise(count = n())

#Count number of PD cases who are missing time to event data
PD_merged %>% 
  filter(is.na(timeToEvent_death)) %>% 
  summarise(count = n())

#Excluding individuals who were only identified as PD from the death records
PD_merged %>% 
  filter(!is.na(timeToEvent_death)) %>% 
  #filter(PD_status == "incident") %>% 
  group_by(event_death) %>% 
  summarise(count = n(),
            mean_time = mean(timeToEvent_death),
            sd_time = sd(timeToEvent_death))

#Grouped by PD status
#Excluding individuals who were only identified as PD from the death records
PD_merged %>% 
  filter(!is.na(timeToEvent_death)) %>% 
  group_by(PD_status, event_death) %>% 
  summarise(count = n(),
            mean_time = mean(timeToEvent_death),
            sd_time = sd(timeToEvent_death))

#Get the earliest and latest PD date of death
PD_merged %>% 
  summarise(min = min(date_of_death, na.rm = TRUE),
            max = max(date_of_death, na.rm = TRUE))

#Mean age at diagnosis excluding incident cases just identified at death
PD_merged %>% 
  filter(!is.na(timeToEvent_death)) %>% 
  group_by(PD_status) %>% 
  summarise(count = n(),
            mean_aad = mean(age_diagnosis, na.rm = TRUE),
            sd_aad = sd(age_diagnosis),
            mean_ageentry = mean(age_entry),
            sd_ageentry = sd(age_entry),
            mean_dd = mean(disease_duration_onset, na.rm = TRUE),
            sd_dd = sd(disease_duration_onset, na.rm = TRUE))

#---Create survival object for mortality/survival---####

surv_object_mortality <- Surv(time = PD_merged$timeToEvent_death, event = PD_merged$event_death)

#---Plot Kaplan-Meier curve for gender vs. mortality---####

fit_gender_mortality <- survfit(surv_object_mortality ~ gender, data = PD_merged)

summary(fit_gender_mortality)

#Plot Kaplan Meier curve for gender vs. mortality
ggsurvplot(fit_gender_mortality, data = PD_merged, pval = TRUE) +
  ggsave("plots/UKB_mortality_gender.png")

#Fit Cox Proportional Hazards model for mortality vs. gender
fit_gender_mortality_coxph = coxph(surv_object_mortality ~ gender, data = PD_merged)
summary(fit_gender_mortality_coxph)

#Fit Cox Proportional Hazards model for mortality vs. gender with covariates
fit_gender_mortality_coxph_covars = coxph(surv_object_mortality ~ gender + age_diagnosis, data = PD_merged)
summary(fit_gender_mortality_coxph_covars)

#---Export clinical data for GWAS survival analysis of mortality---####

#Need to export: event, time to event, age at onset, gender
UKB_export_mortality <- PD_merged %>% 
  select(eid, event_death, timeToEvent_death, age_diagnosis, gender, PD_status, PDdeath) %>% 
  mutate(FID = paste(eid, "_", eid, sep = "")) %>% 
  mutate(IID = FID) %>% 
  select(FID, IID, event_death, timeToEvent_death, age_diagnosis, gender, PD_status, PDdeath) %>% #rearrange order of variables
  filter(!is.na(timeToEvent_death)) #Remove individuals who are missing time to event

#Summarise incident vs. prevalent cases after removing individuals who were only identified at death
UKB_export_mortality %>% 
  group_by(PD_status) %>% 
  summarise(count = n())

#Write as output
filename_mortality <- paste("outputs/UKB_survival_mortality_", today(), ".txt", sep ="")

write.table(UKB_export_mortality, filename_mortality,
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

#Export just list of PD individuals (including the ones just from death data)
UKB_export_PD_IID_FID <- PD_merged %>% 
  select(eid) %>% 
  mutate(FID = eid) %>% 
  rename(IID = eid)

write.table(UKB_export_PD_IID_FID, "outputs/UKB_PDcases_IDs.txt",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")




########## EXPORT CLINICAL DATA FOR MERGE W OTHER COHORTS ##########
#---Export data for merge with other cohorts---####

UKB_export_all <- PD_merged %>%
  filter(PD_status!="undefined") %>%
  mutate(FID = paste(eid, "_", eid, sep = "")) %>% 
  mutate(IID = FID) %>% 
  rename(age_onset_imput = age_diagnosis,
         disease_duration_onset_imput = disease_duration_onset) %>% 
  select(FID, IID, event_death, timeToEvent_death,
         age_onset_imput, gender, PD_status, PDdeath, disease_duration_onset_imput) %>% 
  mutate(cohort = ifelse(PD_status == "incident", "UKB_incident",
                         ifelse(PD_status == "prevalent", "UKB_prevalent", NA))) %>% 
  select(-PD_status) %>%
  mutate(event_HY3 = NA, 
         timeToEvent_HY3 = NA,
         event_dementia = NA, 
         timeToEvent_dementia = NA,
         BL_HY2 = NA, 
         subtype_cat = NA, 
         V1_MOCA_total_adj = NA, 
         years_education_bin = NA)


filename_all <- paste("outputs/UKB_survival_all_", today(), ".txt", sep ="")

write.table(UKB_export_all, filename_all,
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

