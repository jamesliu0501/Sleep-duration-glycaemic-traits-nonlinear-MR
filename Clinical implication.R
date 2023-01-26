# This file presents the R code of the clinical implication
library(data.table)
library(stringi)
library(ff)


# Clinical implication 1
# The estimations were based the effect sizes of the treatment effects in the main self-reported estimates (5 stratum) on HbA1c
# We apply 20 minutes change which is equal to ~ 4 SD change of the genetic variant 
# 1SD = 5.4 and Beta = 0.015 hours per 1 unit change of the allele score, where 1SD = 0.08 hours / ~ 4.8 mins



# Import data
data<- fread(file="data.csv",
             fill=T, header = T)

# subset the data who sleep <= 7 hours/day
data_sleep<-subset(data, data$sleepduration_sr<=7)


# The estimated effect of each stratum 
# mean 5.98hr:  -1.39 (-1.97 to -0.81) mmol/mol / hr (stratum 1, the effect applied to who sleep < 6 hours per day)
# mean 6.70hr:  -0.78 (-1.32 to -0.24) mmol/mol / hr (stratum 2)
# mean 7.17hr:  -0.75 (-1.27 to -0.22) mmol/mol / hr (stratum 3, the effect applied to who sleep >= 6 and <= 7 hours per day)
# mean 7.64hr:  -0.72 (-1.30 to -0.14) mmol/mol / hr (stratum 4)

# before the sleep treatment
data_sleep$hba1c_48<-ifelse(data_sleep$hba1c1>=48, 1, 0)

# among who were <6hr, we lower their HbA1c level by 1.31 * 1/3
data_sleep$hba1c_treated<-ifelse(data_sleep$sleepduration_sr < 6, data_sleep$hba1c1 - 1.31 * 1/3, data_sleep$hba1c1)

# among who were >=6hr but =<7hr (), we lower their HbA1c level by 0.75 * 0.083 mmol/mol
data_sleep$hba1c_treated<-ifelse(data_sleep$sleepduration_sr <= 7 & data_sleep$sleepduration_sr>=6, data_sleep$hba1c1 - 0.75*1/3, data_sleep$hba1c_treated)

# after sleep treatment
data_sleep$hba1c_treated_48<-ifelse(data_sleep$hba1c_treated>=48, 1, 0)

# absolute change of diabetes
nchange<-hba1c_48[2] - hba1c_treated_48[2]

# The conditional percentage change 
cp<- (nchange /  hba1c_48[2])*100









# Implication 2 
# Import data
data<- fread(file="data.csv",
             fill=T, header = T)

# The estimated effect of each stratum 
# mean 5.98hr:  -1.39 (-1.97 to -0.81) mmol/mol / hr (stratum 1, the effect applied to who sleep < 6 hours per day)
# mean 6.70hr:  -0.78 (-1.32 to -0.24) mmol/mol / hr (stratum 2)
# mean 7.17hr:  -0.75 (-1.27 to -0.22) mmol/mol / hr (stratum 3, the effect applied to who sleep >= 6 and <= 7 hours per day)
# mean 7.64hr:  -0.72 (-1.30 to -0.14) mmol/mol / hr (stratum 4)


# Calculate the proportion of diabetes in the UKB before sleep treatment, regardless of sleep duration
data$hba1c_48<-ifelse(data$hba1c1>=48, 1, 0)
hba1c_48<-table(data$hba1c_48)
p0<-hba1c_48[2] / (hba1c_48[1] + hba1c_48[2]) 


# Among who were <6hr, we lower their HbA1c level by 1.31 * 20 / 60 mmol/mol
data$hba1c_treated<-ifelse(data$sleepduration_sr < 6, data$hba1c1 - 1.31 * 1/3, data$hba1c1)

# Among who were >=6hr but =<7hr, we lower their HbA1c level by 0.75 * 20/ 60 mmol/mol
data$hba1c_treated<-ifelse(data$sleepduration_sr <= 7 & data$sleepduration_sr>=6, data$hba1c1 - 0.75*1/3, data$hba1c_treated)


# Calculate the proportion of diabetes in the UKB after sleep treatment
data$hba1c_treated_48<-ifelse(data$hba1c_treated>=48, 1, 0)
hba1c_treated_48<-table(data$hba1c_treated_48)

p1<-hba1c_treated_48[2] / (hba1c_treated_48[1] + hba1c_treated_48[2])

# Calculate the proportion change
preduce<- (p0 - p1)

# We then used a parametric bootstrap to obtain the corresponding 95% confidence interval
SE1        = (-0.81 - -1.39)/1.96*1/3   # standard error 
SE2        = (-0.22 - -0.75)/1.96*1/3   # standard error 

# Vector of HbA1cs for whole sample
## sleepduration < 6hr  I1 = 1, I2 = 0 
## sleepduration = 6, 7hr  I1 = 0, I2 = 1
## sleepduration > 7hr  I1 = 0, I2 = 0
HbA1c     = data$hba1c1          
I1<-ifelse(data$sleepduration_sr<6, 1, 0)
I2<-ifelse(data$sleepduration_sr<=7 & data$sleepduration_sr>=6, 1, 0)


# Calculate the 95% CI
P0        = p0                       
N         = hba1c_treated_48[1] + hba1c_treated_48[2]
P1 = NULL
P01 = NULL
set.seed(12345)
for(i in 1:1000){
  bstar1 =  rnorm(1,-1.39*1/3,SE1)
  bstar2 =  rnorm(1,-0.75*1/3,SE2)
  HbA1cStar = HbA1c + I1*bstar1 + I2*bstar2
  HbA1cStar = na.omit(HbA1cStar)
  P1[i]            = length(HbA1cStar[HbA1cStar>=48])/N
  P01[i]         = P0 - P1[i]
}
SE     = sd(P01)
CI_PI = P01 +c(-1.96,1.96)*SE

P_difference<- P01[1]
CI_PI<-P_difference + c(-1.96,1.96)*SE


# Applying to 37% of the England and Wales population (i.e., 22 million out of 60 million in total) are in 40 - 70 age rang.
n =  22000000 * P_difference
nci = 22000000 * CI_PI


