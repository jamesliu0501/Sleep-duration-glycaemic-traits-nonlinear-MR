# This file includes the R scripts of the main analysis (piecewise linear method (PLM) with doubly-ranked stratification method)
# The details of the nonlinear Mendelian randomization should be referred to https://github.com/amymariemason/SUMnlmr 
# The following scripts are taking the analysis of self-reported (SR) sleep duration with HbA1c (5 strata) as an example. The 3- and 10- strata (and with glucose) cases are similar. 
devtools::install_github("amymariemason/SUMnlmr", force = TRUE)
install.packages("cachem")
library(cachem)
library(SUMnlmr)
library(data.table)
library(ggplot2)
library(ff)



# Import the data 
data<-fread(file="data.csv",  fill=T, header = T)


# Select the related variables in the model
data_sr_hba1c<-select(data, projectid, sleepduration_sr, sleepduration_pl8_unweighted, hba1c1,
                      sex, age_recruitment, assessment_centre1, chip, 
                      PC1,PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                      PC11, PC12, PC13, PC14, PC15, PC16)

# Excluding NA rows
data_sr_hba1c<-na.omit(data_sr_hba1c)

# To creat a dummy variable for a categolorical variable (i.e., centre)
centre <- model.matrix(~as.factor(assessment_centre1), data = data_sr_hba1c)
centre <- centre[,-1]

# A list of the covariates
cov<-cbind(data_sr_hba1c$age_recruitment, centre, data_sr_hba1c$sex, as.factor(data_sr_hba1c$chip),
           data_sr_hba1c$PC1,data_sr_hba1c$PC2,data_sr_hba1c$PC3,data_sr_hba1c$PC4,data_sr_hba1c$PC5,data_sr_hba1c$PC6,data_sr_hba1c$PC7,data_sr_hba1c$PC8,data_sr_hba1c$PC9,data_sr_hba1c$PC10,
           data_sr_hba1c$PC11,data_sr_hba1c$PC12,data_sr_hba1c$PC13,data_sr_hba1c$PC14,data_sr_hba1c$PC15,data_sr_hba1c$PC16)


# Main analysis (doubly-ranked stratification)
# Generate the summary data into 5 strata with the doubly-reanked stratification method (the cases of 3 and 10 strata are similiar)
set.seed(12345) 
summ_sr_hba1c_5s<-create_nlmr_summary(y = data_sr_hba1c$hba1c1,
                                      x = data_sr_hba1c$sleepduration_sr,
                                      g = data_sr_hba1c$sleepduration_pl8_unweighted,
                                      covar =cov,
                                      family = "gaussian",
                                      controlsonly = FALSE,
                                      strata_method = "ranked" ,
                                      q = 5)


# The piecewise linear model (set 7 hours as the reference group)
PLM_sr_hba1c_5s <-with(summ_sr_hba1c_5s$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                                   ci="bootstrap_se",
                                                                   nboot=1000, 
                                                                   fig=TRUE,
                                                                   family="gaussian",
                                                                   ref = 7,
                                                                   ci_fig="ribbon"))



# Present the figure 
plot_PLM_sr_hba1c_5s <- PLM_sr_hba1c_5s$figure + ggtitle("SR sleep duration vs HbA1c (PLM 5s_rank)") + theme(plot.title = element_text(size=25))

# Summary of the estimates 
sum_PLM_sr_hba1c_5s<-summary(PLM_sr_hba1c_5s)





# Sensitivity analysis (residual stratification)
# Generate the summary data by iv-free exposure into 3 strata (for self-reported sleep duration, only 3 strata is stratified; for accelerometer-derived sleep duration, 3, 5, and 10 strata are stratified)
summ_sr_hba1c_3s<-create_nlmr_summary(y = data_sr$hba1c1,
                                      x = data_sr$sleepduration_sr,
                                      g = data_sr$sleepduration_pl8_unweighted,
                                      covar =cov,
                                      family = "gaussian",
                                      controlsonly = FALSE,
                                      q = 3,
                                      strata_method = "residual" ,
                                      extra_statistics = TRUE)                                      

# Piecewise linear model (set 7 hours as the reference)
PLM_sr_hba1c_3s <-with(summ_sr_hba1c_3s$summary, piecewise_summ_mr(by, bx, byse, bxse, xmean, xmin,xmax, 
                                                                   ci="bootstrap_se",
                                                                   nboot=1000, 
                                                                   fig=TRUE,
                                                                   family="gaussian",
                                                                   ref = 7,
                                                                   ci_fig="ribbon"))

# Present the figure 
plot_PLM_sr_hba1c_3s <- PLM_sr_hba1c_3s$figure + ggtitle("SR sleep duration vs HbA1c (PLM 3s-residual)")

# Summary of the estimates 
sum_PLM_sr_hba1c_3s<-summary(PLM_sr_hba1c_3s)

