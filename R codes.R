rm(list=ls())
library(readxl)
library(dilutBMS2)
library(BMS)
library(fixest)
library(corrplot)
library(metafor)

library(multiwayvcov)
library(foreign)
library(multcomp)
library(ggplot2)
library(dplyr)
library(forcats)
library(AER)
library(puniform)
library(Hmisc)
library(DescTools)
library(plm)
library(fwildclusterboot)
library(fixest)
library("RoBMA")
library("data.table")
library(MetaStudies)
library(phack)
library(tidyverse)


Data <- read_excel("C:/Users/USER PC/Documents/Data.xlsx", sheet = "All")
View(Data)




###Overnight interbank reference rate sample

########################################################
########### WAAP by  Ioannidis et al. (2017) ###########
########################################################

#Upward pass-through
Hike <- Data[Data$hike == 1 & Data$inter_overnight ==1, ]
WAAP_true_ef<-sum(Hike$beta)/nrow(Hike)
WAAP_data<-Hike[WAAP_true_ef/2.8>Hike$se,]

#The weights
Inv_SE = 1/(WAAP_data$se*WAAP_data$se)
WAAP_data$Inv_SE <- Inv_SE

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster="study_id", weights = Inv_SE, data = WAAP_data)    
summary(WAAP)



#Downward pass-through
Cut <- Data[Data$cut == 1 & Data$inter_overnight ==1, ]
WAAP_true_ef<-sum(Cut$beta)/nrow(Cut)
WAAP_data<-Cut[WAAP_true_ef/2.8>Cut$se,]

#The weights
Inv_SE = 1/(WAAP_data$se*WAAP_data$se)
WAAP_data$Inv_SE <- Inv_SE

WAAP<-feols(beta ~ se, ssc=ssc(adj = FALSE, cluster.adj=FALSE), cluster="study_id", weights = Inv_SE, data = WAAP_data)    
summary(WAAP)






#################################################
#### Stem based technique by Furukawa (2021) ####
#################################################

#visit https://github.com/Chishio318/stem-based_method and download the R script (stem_method.R)
#paste the script in your working directory and run the file


#Upward pass-through
Hike <- Data[Data$hike == 1 & Data$inter_overnight ==1, ]
source("stem_method.R")
median_data <- data_median(Hike, "study_id", "beta", "se")
stem_results <- stem(median_data$coefficient, median_data$standard_error, param)
View(stem_results$estimates)


#Downward pass-through
Cut <- Data[Data$cut == 1 & Data$inter_overnight ==1, ]
source("stem_method.R")
median_data <- data_median(Cut, "study_id", "beta", "se")
stem_results <- stem(median_data$coefficient, median_data$standard_error, param)
View(stem_results$estimates)





###################################################################################
#### Endogenous Kink model by Bom and Rachinger (2019) using the Meta-Analysis ####
#### Instrumental Variables Estimator package by Havrankova et al. (2023)      ####
###################################################################################

#Download the MAIVE package from: http://meta-analysis.cz/maive/
#Copy the relevant subset of the data and paste into the inputdata excel file in the folder "R code for MAIVE", 
#Open the file maive.R in R studio:
#Select the following parameters:

#method <- 4  #choose the endogenous kink model
#weight <- 1 #normal weights
#instrument <- 0 #no instrument; sample size is a weak instrument in this sample
#studylevel <- 2 #study-level clustered standard errors

#Select (ctrl + a) and run the code
#Repeat the procedure for remaining subsets of the dataset





###############################################################
########### P-hacking tests by Elliot et al. (2021) ###########
###############################################################


#Upward pass-through

#Create backup
Backup <- Data

Data <- Backup[Backup$hike == 1 & Backup$inter_overnight ==1, ]

# Derounding:
Data$added <- 0.005
Data$Ones <- 1

p_upper <- Data$p + Data$added
p_lower <- Data$p - Data$added

p.deround = runif(Data$Ones, p_lower, p_upper)

Data$p.deround <- p.deround

study_id <- Data$study_id

# p-hacking tests
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.05, p_max=0.2, J=5, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.05, p_max=0.2, J=5)

Obs <- Data[Data$p.deround>=0.05 & Data$p.deround<=0.2,]
count(Obs)

#Return backed-up data
Data <- Backup



#Downward pass-through

#Create backup
Backup <- Data

Data <- Backup[Backup$cut == 1 & Backup$inter_overnight ==1, ]

# Derounding:
Data$added <- 0.005
Data$Ones <- 1

p_upper <- Data$p + Data$added
p_lower <- Data$p - Data$added

p.deround = runif(Data$Ones, p_lower, p_upper)

Data$p.deround <- p.deround

study_id <- Data$study_id

# p-hacking tests
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.05, p_max=0.2, J=5, K=2, use_bound=TRUE)
phack_test_cox_shi(Data$p.deround, Data$study_id, p_min=0.05, p_max=0.2, J=5)

Obs <- Data[Data$p.deround>=0.05 & Data$p.deround<=0.2,]
count(Obs)

#Return backed-up data
Data <- Backup








##clear environment before running BMA models
rm(list=ls())

## Bayesian model averaging (Overnight interbank reference rate sample)

# Importing data containing overnight interbank reference rates only

Data <- read_excel("C:/Users/USER PC/Documents/Data.xlsx", sheet = "overnight_interbank")
View(Data)


#Extracting the interacted variables

Data$CBI_Mpol <- Data$CBI*Data$M_policy
Data$CBI_Dev <- Data$CBI*Data$Dev #not used
Data$CBI_Ex <- Data$CBI*Data$Ex_regime#not used
Data$CBI_3 <- Data$CBI*Data$M_policy*Data$Dev#not used
Data$CBI_4 <- Data$CBI*Data$Dev#not used





#Weighted by inverse number of estimates

#re-extracting the variables
IRPT <- Data[, c("beta", "se", "hike", "cut", "Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "PreCrises", "PostCrises_PreCovid", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "Mcred")]

#Construct Weights:
Data$Weight_est <- 1/Data$est_no

#Weight the variables:
IRPT <- IRPT*Data$Weight_est

#running the model
bma_irpt2 = dilutBMS2::bms(IRPT, burn=1000000, iter=2000000, mprior="dilut", g="UIP")






#Determinants of asymmetry

#Weighted by inverse number of estimates

#re-extracting the variables

IRPT1 <- Data[, c("beta", "se", "hike", "cut","Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite")]

IRPT2 <- Data[, c("PreCrises", "PostCrises_PreCovid", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "Mcred")]


IRPT2 <- IRPT2*Data$hike

IRPT <- cbind(IRPT1, IRPT2)


#Construct Weights:
Data$Weight_est <- 1/Data$est_no

#Weight the variables:
IRPT <- IRPT*Data$Weight_est

#running the model
bma_irpt7 = dilutBMS2::bms(IRPT, burn=1000000, iter=2000000, mprior="dilut", g="UIP")





##############################
###### Implied estimates #####
##############################


Data$Inv_est <- 1/(Data$est_no*Data$est_no)
Data <- na.omit(Data)


##############################
###### Overall estimates #####
##############################

#Overall: Upward long run
new.dat <- data.frame(
  se=0, 
  hike =1, 
  cut=0, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=mean(Data$Firms), 
  Consumer_Household=mean(Data$Consumer_Household), 
  Mortgage=mean(Data$Mortgage), 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')


#Overall: Downward long run
new.dat <- data.frame(
  se=0, 
  hike =0, 
  cut=1, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=mean(Data$Firms), 
  Consumer_Household=mean(Data$Consumer_Household), 
  Mortgage=mean(Data$Mortgage), 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')



##############################
## Business loan estimates ###
##############################

#Overall: Upward long run
new.dat <- data.frame(
  se=0, 
  hike =1, 
  cut=0, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=1, 
  Consumer_Household=0, 
  Mortgage=0, 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')


#Overall: Downward long run
new.dat <- data.frame(
  se=0, 
  hike =0, 
  cut=1, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=1, 
  Consumer_Household=0, 
  Mortgage=0, 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')




##############################
## Consumer loan estimates ###
##############################

#Overall: Upward long run
new.dat <- data.frame(
  se=0, 
  hike =1, 
  cut=0, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=0, 
  Consumer_Household=1, 
  Mortgage=0, 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')


#Overall: Downward long run
new.dat <- data.frame(
  se=0, 
  hike =0, 
  cut=1, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=0, 
  Consumer_Household=1, 
  Mortgage=0, 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')




##############################
## Mortgage loan estimates ###
##############################

#Overall: Upward long run
new.dat <- data.frame(
  se=0, 
  hike =1, 
  cut=0, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=0, 
  Consumer_Household=0, 
  Mortgage=1, 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')


#Overall: Downward long run
new.dat <- data.frame(
  se=0, 
  hike =0, 
  cut=1, 
  Annual=0, 
  Quarterly=0, 
  Monthly=1, 
  IMF=mean(Data$IMF), 
  ECB=mean(Data$ECB), 
  Time_Series=0, 
  Microdata=1, 
  PreCrises=mean(Data$PreCrises), 
  PostCrises_PreCovid=mean(Data$PostCrises_PreCovid), 
  Short_run=0, 
  Lag_length=0, 
  Infl_con=0, 
  OLS=mean(Data$OLS), 
  Peer=1, 
  IMPF=max(Data$IMPF), 
  A_Cite=max(Data$A_Cite), 
  LP_year=max(Data$LP_year),
  Central_banker=mean(Data$Central_banker), 
  Firms=0, 
  Consumer_Household=0, 
  Mortgage=1, 
  S_loan=mean(Data$S_loan), 
  Imp_Open=mean(Data$Imp_Open), 
  FDI_Open=mean(Data$FDI_Open), 
  Dev=mean(Data$Dev), 
  Cap=mean(Data$Cap), 
  Stock_turn=mean(Data$Stock_turn), 
  CBI_Mpol=mean(Data$CBI_Mpol), 
  Ex_regime=mean(Data$Ex_regime), 
  Eurozone=mean(Data$Eurozone), 
  Net_save=mean(Data$Net_save),
  Mcred=mean(Data$Mcred))

predict(bma_irpt2, newdata = new.dat, interval = 'confidence')


#OLS confidence intervals
OLSCONF <- feols(beta ~ se + hike + cut + Annual + Quarterly + Monthly + IMF + ECB + Time_Series + Microdata + PreCrises + PostCrises_PreCovid + Short_run + Lag_length + Infl_con + OLS + Peer + IMPF + A_Cite + LP_year + Central_banker + Firms + Consumer_Household + Mortgage + S_loan + Imp_Open + FDI_Open + Dev + Cap + Stock_turn + CBI_Mpol + Ex_regime + Eurozone + Net_save + Mcred,  ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Inv_est, cluster = "study_id", panel.id = "study_id",  data = Data)
predict(OLSCONF, newdata = new.dat, interval = 'confidence')



#Diagnostics and other output for Baseline BMA model

#Renaming the variables

#Baseline BMA model output:
bma_irpt2$reg.names[[1]] <- "Standard error"
bma_irpt2$reg.names[[2]] <- "Upward pass-through"
bma_irpt2$reg.names[[3]] <- "Downward pass-through"
bma_irpt2$reg.names[[4]] <- "Annual data"
bma_irpt2$reg.names[[5]] <- "Quarterly data"
bma_irpt2$reg.names[[6]] <- "Monthly data"
bma_irpt2$reg.names[[7]] <- "IMF data"
bma_irpt2$reg.names[[8]] <- "ECB data"
bma_irpt2$reg.names[[9]] <- "Time series"
bma_irpt2$reg.names[[10]] <- "Micro-level data"
bma_irpt2$reg.names[[11]] <- "Pre-GFC data"
bma_irpt2$reg.names[[12]] <- "Post-GFC/Pre-COVID data"
bma_irpt2$reg.names[[13]] <- "Short-run"
bma_irpt2$reg.names[[14]] <- "Lag length"
bma_irpt2$reg.names[[15]] <- "Simultaneity bias"
bma_irpt2$reg.names[[16]] <- "OLS"
bma_irpt2$reg.names[[17]] <- "Publication year"
bma_irpt2$reg.names[[18]] <- "Central bank author"
bma_irpt2$reg.names[[19]] <- "Peer-reviewed journal"
bma_irpt2$reg.names[[20]] <- "Impact factor"
bma_irpt2$reg.names[[21]] <- "Annual citations"
bma_irpt2$reg.names[[22]] <- "Business loan rate"
bma_irpt2$reg.names[[23]] <- "Consumer loan rate"
bma_irpt2$reg.names[[24]] <- "Mortgage loan rate"
bma_irpt2$reg.names[[25]] <- "Short-term loan rate"
bma_irpt2$reg.names[[26]] <- "Trade openness"
bma_irpt2$reg.names[[27]] <- "FDI Openness"
bma_irpt2$reg.names[[28]] <- "Development status"
bma_irpt2$reg.names[[29]] <- "Market capitalization"
bma_irpt2$reg.names[[30]] <- "Stock turnover ratio"
bma_irpt2$reg.names[[31]] <- "Central bank independence (IT)"
bma_irpt2$reg.names[[32]] <- "Floating exchange rate"
bma_irpt2$reg.names[[33]] <- "Euro area member"
bma_irpt2$reg.names[[34]] <- "Net national savings"
bma_irpt2$reg.names[[35]] <- "Credit to private sector"



#Baseline BMA model inclusion graph
image(bma_irpt2,order.by.pip = TRUE, cex.axis = 0.9, main=NULL, xlab = "", do.par = TRUE)


#Baseline BMA model size and convergence
plot(bma_irpt2, include.legend=TRUE)


#Baseline BMA model fitted intercept
coef(bma_irpt2, include.constant = TRUE)


#Baseline BMA model summary
summary(bma_irpt2)


#Coefficient distributions: plotted sequentially
density(bma_irpt2, reg = 2, main = "PIP: 63.64%") #Upward pass-through
density(bma_irpt2, reg = 26, main = "PIP: 91.16%") #Trade (import) openness
density(bma_irpt2, reg = 34, main = "PIP: 100%") #Net national savings
density(bma_irpt2, reg = 35, main = "PIP: 80.87%") #Credit to private sector




#Correlation matrix
bma_COR <- Data[, c("se", "hike", "cut", "Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "PreCrises", "PostCrises_PreCovid", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Mcred", "Cap", "Stock_turn", "CBI_Mpol", "Net_save", "Ex_regime", "Eurozone")]

names(bma_COR)[names(bma_COR) =="se"] <- "Standard error"
names(bma_COR)[names(bma_COR) =="hike"] <- "Upward pass-through"
names(bma_COR)[names(bma_COR) =="cut"] <- "Downward pass-through"
names(bma_COR)[names(bma_COR) =="Annual"] <- "Annual data"
names(bma_COR)[names(bma_COR) =="Quarterly"] <- "Quarterly data"
names(bma_COR)[names(bma_COR) =="Monthly"] <- "Monthly data"
names(bma_COR)[names(bma_COR) =="IMF"] <- "IMF data"
names(bma_COR)[names(bma_COR) =="ECB"] <- "ECB data"
names(bma_COR)[names(bma_COR) =="Time_Series"] <- "Time series"
names(bma_COR)[names(bma_COR) =="Microdata"] <- "Micro-level data"
names(bma_COR)[names(bma_COR) =="PreCrises"] <- "Pre-GFC"
names(bma_COR)[names(bma_COR) =="PostCrises_PreCovid"] <- "Post-GFC/Pre-COVID"
names(bma_COR)[names(bma_COR) =="Short_run"] <- "Short-run"
names(bma_COR)[names(bma_COR) =="Lag_length"] <- "Lag length"
names(bma_COR)[names(bma_COR) =="Infl_con"] <- "Simultaneity bias"
names(bma_COR)[names(bma_COR) =="OLS"] <- "OLS"
names(bma_COR)[names(bma_COR) =="LP_year"] <- "Publication year"
names(bma_COR)[names(bma_COR) =="Central_banker"] <- "Central bank author"
names(bma_COR)[names(bma_COR) =="Peer"] <- "Peer-reviewed journal"
names(bma_COR)[names(bma_COR) =="IMPF"] <- "Impact factor"
names(bma_COR)[names(bma_COR) =="A_Cite"] <- "Annual citations"
names(bma_COR)[names(bma_COR) =="Firms"] <- "Business loan"
names(bma_COR)[names(bma_COR) =="Consumer_Household"] <- "Consumer loan"
names(bma_COR)[names(bma_COR) =="Mortgage"] <- "Mortgage loan"
names(bma_COR)[names(bma_COR) =="S_loan"] <- "Short-term loan"
names(bma_COR)[names(bma_COR) =="Imp_Open"] <- "Trade openness"
names(bma_COR)[names(bma_COR) =="FDI_Open"] <- "FDI openness"
names(bma_COR)[names(bma_COR) =="Dev"] <- "Development status"
names(bma_COR)[names(bma_COR) =="Cap"] <- "Market capitalization"
names(bma_COR)[names(bma_COR) =="Stock_turn"] <- "Stock turnover ratio"
names(bma_COR)[names(bma_COR) =="CBI_Mpol"] <- "Central bank independence (IT)"
names(bma_COR)[names(bma_COR) =="Ex_regime"] <- "Floating exchange rate"
names(bma_COR)[names(bma_COR) =="Eurozone"] <- "Euro area member"
names(bma_COR)[names(bma_COR) =="Net_save"] <- "Net national savings"
names(bma_COR)[names(bma_COR) =="Mcred"] <- "Credit to private sector"


bma_COR <- na.omit(bma_COR)
BMA_COR<-cor(bma_COR)
corrplot(BMA_COR, type="lower", tl.cex = 0.55, number.cex = 0.2, sig.level = 0.01, tl.col="black", addCoef.col = "black", method="color", tl.srt=45)







#Diagnostics and other output for second BMA model (determinants of asymmetry)

#Renaming the variables

#Baseline BMA model output:
bma_irpt7$reg.names[[1]] <- "Standard error"
bma_irpt7$reg.names[[2]] <- "Upward pass-through"
bma_irpt7$reg.names[[3]] <- "Downward pass-through"
bma_irpt7$reg.names[[4]] <- "Annual data"
bma_irpt7$reg.names[[5]] <- "Quarterly data"
bma_irpt7$reg.names[[6]] <- "Monthly data"
bma_irpt7$reg.names[[7]] <- "IMF data"
bma_irpt7$reg.names[[8]] <- "ECB data"
bma_irpt7$reg.names[[9]] <- "Time series"
bma_irpt7$reg.names[[10]] <- "Micro-level data"
bma_irpt7$reg.names[[11]] <- "Pre-GFC data"
bma_irpt7$reg.names[[12]] <- "Post-GFC/Pre-COVID data"
bma_irpt7$reg.names[[13]] <- "Short-run"
bma_irpt7$reg.names[[14]] <- "Lag length"
bma_irpt7$reg.names[[15]] <- "Simultaneity bias"
bma_irpt7$reg.names[[16]] <- "OLS"
bma_irpt7$reg.names[[17]] <- "Publication year"
bma_irpt7$reg.names[[18]] <- "Central bank author"
bma_irpt7$reg.names[[19]] <- "Peer-reviewd journal"
bma_irpt7$reg.names[[20]] <- "Impact factor"
bma_irpt7$reg.names[[21]] <- "Annual citations"
bma_irpt7$reg.names[[22]] <- "Business loan rate"
bma_irpt7$reg.names[[23]] <- "Consumer loan rate"
bma_irpt7$reg.names[[24]] <- "Mortgage loan rate"
bma_irpt7$reg.names[[25]] <- "Short-term loan"
bma_irpt7$reg.names[[26]] <- "Trade openness"
bma_irpt7$reg.names[[27]] <- "FDI Openness"
bma_irpt7$reg.names[[28]] <- "Development status"
bma_irpt7$reg.names[[29]] <- "Market capitalization"
bma_irpt7$reg.names[[30]] <- "Stock turnover ratio"
bma_irpt7$reg.names[[31]] <- "Central bank independence"
bma_irpt7$reg.names[[32]] <- "Floating exchange rate"
bma_irpt7$reg.names[[33]] <- "Euro area member"
bma_irpt7$reg.names[[34]] <- "Net national savings"
bma_irpt7$reg.names[[35]] <- "Credit to private sector"



#Baseline BMA model inclusion graph
image(bma_irpt7,order.by.pip = TRUE, cex.axis = 0.9, main=NULL, xlab = "", do.par = TRUE)


#Baseline BMA model size and convergence
plot(bma_irpt7, include.legend=TRUE)


#Baseline BMA model fitted intercept
coef(bma_irpt7, include.constant = TRUE)


#Baseline BMA model summary
summary(bma_irpt7)




#================================
# Frequentist check 1
#================================

#Frequentist check
IRPT <- Data[, c("beta", "se", "hike", "cut", "Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "PreCrises", "PostCrises_PreCovid", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "Mcred")]

#Construct Weights:
IRPT$Weight_est <- 1/(Data$est_no*Data$est_no)
IRPT$study_id <- Data$study_id
IRPT <- na.omit(IRPT)

freq_ols <-feols(beta ~ se + OLS + Central_banker + Net_save + Annual + Microdata + PostCrises_PreCovid + Infl_con + Monthly + A_Cite + ECB + Imp_Open + Consumer_Household + Mcred + Lag_length + IMPF + hike, ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, data = IRPT)

summary(freq_ols)
coeftest(freq_ols, vcov=vcovHC(freq_ols,type="HC0", cluster = c(IRPT$study_id)))



#================================
# Frequentist check 2
#================================

#Frequentist check

IRPT1 <- Data[, c("beta", "se", "hike", "cut","Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite")]
IRPT2 <- Data[, c("PreCrises", "PostCrises_PreCovid", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "Mcred")]


IRPT2 <- IRPT2*Data$hike

IRPT <- cbind(IRPT1, IRPT2)

#Construct Weights:
IRPT$Weight_est <- 1/(Data$est_no*Data$est_no)
IRPT$study_id <- Data$study_id
IRPT <- na.omit(IRPT)

freq_ols <-feols(beta ~ se + Annual + PreCrises + Imp_Open + Eurozone + Net_save + CBI_Mpol + Cap + OLS + Lag_length + Stock_turn + hike + Infl_con + ECB, ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights=~Weight_est, data = IRPT)

summary(freq_ols)
coeftest(freq_ols, vcov=vcovHC(freq_ols,type="HC0", cluster = c(IRPT$study_id)))






##clear environment again before running BMA models
rm(list=ls())

## Extended BMA model (full sample)

# Re-importing data containing all reference rates

Data <- read_excel("C:/Users/USER PC/Documents/Data.xlsx", sheet = "All")
View(Data)

#Extracting the interacted variables

Data$CBI_Mpol <- Data$CBI*Data$M_policy
Data$CBI_Dev <- Data$CBI*Data$Dev #not used
Data$CBI_Ex <- Data$CBI*Data$Ex_regime#not used
Data$CBI_3 <- Data$CBI*Data$M_policy*Data$Dev#not used
Data$CBI_4 <- Data$CBI*Data$Dev#not used





#Weighted by inverse number of estimates

#re-extracting the variables
IRPT <- Data[, c("beta", "se", "hike", "cut", "Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "PreCrises", "PostCrises_PreCovid", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "Mcred", "inter_overnight", "inter_3months", "CB_Discount", "CB_overnight")]

#Construct Weights:
Data$Weight_est <- 1/Data$est_no

#Weight the variables:
IRPT <- IRPT*Data$Weight_est

#running the model
bma_irpt8 = dilutBMS2::bms(IRPT, burn=1000000, iter=2000000, mprior="dilut", g="UIP")






###Diagnostics and other output for Baseline BMA model

#Renaming the variables

bma_irpt8$reg.names[[1]] <- "Standard error"
bma_irpt8$reg.names[[2]] <- "Upward pass-through"
bma_irpt8$reg.names[[3]] <- "Downward pass-through"
bma_irpt8$reg.names[[4]] <- "Annual data"
bma_irpt8$reg.names[[5]] <- "Quarterly data"
bma_irpt8$reg.names[[6]] <- "Monthly data"
bma_irpt8$reg.names[[7]] <- "IMF data"
bma_irpt8$reg.names[[8]] <- "ECB data"
bma_irpt8$reg.names[[9]] <- "Time series"
bma_irpt8$reg.names[[10]] <- "Micro-level data"
bma_irpt8$reg.names[[11]] <- "Pre-GFC data"
bma_irpt8$reg.names[[12]] <- "Post-GFC/Pre-COVID data"
bma_irpt8$reg.names[[13]] <- "Short-run"
bma_irpt8$reg.names[[14]] <- "Lag length"
bma_irpt8$reg.names[[15]] <- "Simultaneity bias"
bma_irpt8$reg.names[[16]] <- "OLS"
bma_irpt8$reg.names[[17]] <- "Publication year"
bma_irpt8$reg.names[[18]] <- "Central bank author"
bma_irpt8$reg.names[[19]] <- "Peer-reviewed journal"
bma_irpt8$reg.names[[20]] <- "Impact factor"
bma_irpt8$reg.names[[21]] <- "Annual citations"
bma_irpt8$reg.names[[22]] <- "Business loan rate"
bma_irpt8$reg.names[[23]] <- "Consumer loan rate"
bma_irpt8$reg.names[[24]] <- "Mortgage loan rate"
bma_irpt8$reg.names[[25]] <- "Short-term loan rate"
bma_irpt8$reg.names[[26]] <- "Trade openness"
bma_irpt8$reg.names[[27]] <- "FDI Openness"
bma_irpt8$reg.names[[28]] <- "Development status"
bma_irpt8$reg.names[[29]] <- "Market capitalization"
bma_irpt8$reg.names[[30]] <- "Stock turnover ratio"
bma_irpt8$reg.names[[31]] <- "Central bank independence (IT)"
bma_irpt8$reg.names[[32]] <- "Floating exchange rate"
bma_irpt8$reg.names[[33]] <- "Euro area member"
bma_irpt8$reg.names[[34]] <- "Net national savings"
bma_irpt8$reg.names[[35]] <- "Non-interest banking"
bma_irpt8$reg.names[[36]] <- "Credit to private sector"
bma_irpt8$reg.names[[37]] <- "Interbank overnight"
bma_irpt8$reg.names[[38]] <- "Interbank 3 months"
bma_irpt8$reg.names[[39]] <- "CB discount rate"
bma_irpt8$reg.names[[40]] <- "CB overnight rate"




#Baseline BMA model inclusion graph
image(bma_irpt8,order.by.pip = TRUE, cex.axis = 0.8, main=NULL, xlab = "", do.par = TRUE)


#Baseline BMA model size and convergence
plot(bma_irpt8, include.legend=TRUE)


#Baseline BMA model fitted intercept
coef(bma_irpt8, include.constant = TRUE)


#Baseline BMA model summary
summary(bma_irpt8)





#Correlation matrix

bma_COR2 <- Data[, c("se", "hike", "cut", "Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "PreCrises", "PostCrises_PreCovid", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "NIB", "Mcred", "inter_overnight", "inter_3months", "CB_Discount", "CB_overnight")]

names(bma_COR2)[names(bma_COR2) =="se"] <- "Standard error"
names(bma_COR2)[names(bma_COR2) =="hike"] <- "Upward pass-through"
names(bma_COR2)[names(bma_COR2) =="cut"] <- "Downward pass-through"
names(bma_COR2)[names(bma_COR2) =="Annual"] <- "Annual data"
names(bma_COR2)[names(bma_COR2) =="Quarterly"] <- "Quarterly data"
names(bma_COR2)[names(bma_COR2) =="Monthly"] <- "Monthly data"
names(bma_COR2)[names(bma_COR2) =="IMF"] <- "IMF data"
names(bma_COR2)[names(bma_COR2) =="ECB"] <- "ECB data"
names(bma_COR2)[names(bma_COR2) =="Time_Series"] <- "Time series"
names(bma_COR2)[names(bma_COR2) =="Microdata"] <- "Micro-level data"
names(bma_COR2)[names(bma_COR2) =="PreCrises"] <- "Pre-GFC"
names(bma_COR2)[names(bma_COR2) =="PostCrises_PreCovid"] <- "Post-GFC/Pre-COVID"
names(bma_COR2)[names(bma_COR2) =="Short_run"] <- "Short-run"
names(bma_COR2)[names(bma_COR2) =="Lag_length"] <- "Lag length"
names(bma_COR2)[names(bma_COR2) =="Infl_con"] <- "Simultaneity bias"
names(bma_COR2)[names(bma_COR2) =="OLS"] <- "OLS"
names(bma_COR2)[names(bma_COR2) =="LP_year"] <- "Publication year"
names(bma_COR2)[names(bma_COR2) =="Central_banker"] <- "Central bank author"
names(bma_COR2)[names(bma_COR2) =="Peer"] <- "Peer-reviewed journal"
names(bma_COR2)[names(bma_COR2) =="IMPF"] <- "Impact factor"
names(bma_COR2)[names(bma_COR2) =="A_Cite"] <- "Annual citations"
names(bma_COR2)[names(bma_COR2) =="Firms"] <- "Business loan"
names(bma_COR2)[names(bma_COR2) =="Consumer_Household"] <- "Consumer loan"
names(bma_COR2)[names(bma_COR2) =="Mortgage"] <- "Mortgage loan"
names(bma_COR2)[names(bma_COR2) =="S_loan"] <- "Short-term loan"
names(bma_COR2)[names(bma_COR2) =="Imp_Open"] <- "Trade openness"
names(bma_COR2)[names(bma_COR2) =="FDI_Open"] <- "FDI openness"
names(bma_COR2)[names(bma_COR2) =="Dev"] <- "Development status"
names(bma_COR2)[names(bma_COR2) =="Cap"] <- "Market capitalization"
names(bma_COR2)[names(bma_COR2) =="Stock_turn"] <- "Stock turnover ratio"
names(bma_COR2)[names(bma_COR2) =="CBI_Mpol"] <- "Central bank independence (IT)"
names(bma_COR2)[names(bma_COR2) =="Ex_regime"] <- "Floating exchange rate"
names(bma_COR2)[names(bma_COR2) =="Eurozone"] <- "Euro area member"
names(bma_COR2)[names(bma_COR2) =="Net_save"] <- "Net national savings"
names(bma_COR2)[names(bma_COR2) =="NIB"] <- "Non-interest banking"
names(bma_COR2)[names(bma_COR2) =="Mcred"] <- "Credit to private sector"
names(bma_COR2)[names(bma_COR2) =="inter_overnight"] <- "Interbank overnight"
names(bma_COR2)[names(bma_COR2) =="inter_3months"] <- "Interbank 3 months"
names(bma_COR2)[names(bma_COR2) =="CB_Discount"] <- "CB discount"
names(bma_COR2)[names(bma_COR2) =="CB_overnight"] <- "CB overnight"



bma_COR2 <- na.omit(bma_COR2)
bma_COR2<-cor(bma_COR2)
corrplot(bma_COR2, type="lower", tl.cex = 0.55, number.cex = 0.2, sig.level = 0.01, tl.col="black", addCoef.col = "black", method="color", tl.srt=45)




#================================
# Frequentist check 3
#================================

#Frequentist check
IRPT <- Data[, c("beta", "se", "hike", "cut", "Annual", "Quarterly", "Monthly", "IMF", "ECB", "Time_Series", "Microdata", "PreCrises", "PostCrises_PreCovid", "Short_run", "Lag_length", "Infl_con", "OLS", "LP_year", "Central_banker", "Peer", "IMPF", "A_Cite", "Firms", "Consumer_Household", "Mortgage", "S_loan", "Imp_Open", "FDI_Open", "Dev", "Cap", "Stock_turn", "CBI_Mpol", "Ex_regime", "Eurozone", "Net_save", "Mcred", "inter_overnight", "inter_3months", "CB_Discount", "CB_overnight")]

#Construct Weights:
IRPT$Weight_est <- 1/(Data$est_no*Data$est_no)
IRPT$study_id <- Data$study_id
IRPT <- na.omit(IRPT)

freq_ols <-feols(beta ~ se + cut + Short_run + OLS + Peer + Net_save + Central_banker + Monthly + Consumer_Household + Mcred + inter_overnight + Annual + Lag_length, ssc=ssc(adj = FALSE, cluster.adj=FALSE), weights = ~Weight_est, data = IRPT)

summary(freq_ols)
coeftest(freq_ols, vcov=vcovHC(freq_ols,type="HC0", cluster = c(IRPT$study_id)))
