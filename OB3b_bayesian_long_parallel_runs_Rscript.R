## R script for:
##  Transgenerational effects of maternal obesogenic diet on appetite: a meta-analysis  
  
#  Step 4 - Bayesian models - setting parallel runs for MCMCglmm models on ML work computer - do in Terminal only!


# IMPORTANT
# Before running this script in the Terminal, Comment/Uncomment sections, as needed
# To run this script paste in the Terminal: source("OB3b_bayesian_long_parallel_runs_Rscript.R")

options(scipen=100)
rm(list=ls())

#getRversion()
install.packages(c("MCMCglmm", "foreach", "doMC"))
library(Matrix)
library(lme4)
library(MCMCglmm)
library(foreach)
library(doMC)
registerDoMC()

# set working directory, as neede: setwd("")

load("data_NC.RData") 
load("data.RData") 
load("VC_matrices_NC.RData") 
load("VC_matrices.RData") 

#prior1 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002))) #inverse-Gamma piror (slightly informative) - the default is non-informative prior
#prior2 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000))) # parameter expanded prior (close to non-informative)
prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
prior3s <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
#prior4 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,fix=1))) #parameter expanded prior with AnivG


### MCMCglmm

#Use strains not species as random effect, because there are 8 strains but only 2 species .
#We will run a separate model with strains as fixed effects.


## Allometrically adjusted food intake - RelIntake2 (RI2)
#-----------------------------------------------------
  
#  Define full dataset variables.
d <- data$RelIntake2_Hd
vofd <- data$RelIntake2_Hd_Var_adj
m <- RelIntake2_matrix
AinvG <- solve(m)
AnivG <- as(AinvG,"dgCMatrix")

# RI2 - Null model for the full data set.

##test runs - without VCV matrix
#m0.0<-MCMCglmm(d~1,random=~Study_ID+Strain,family="gaussian",data=data,verbose=T,mev=vofd,nitt=200000,thin=25,burnin=25000,pr=T) #use stronger prior!
#prior1 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002))) #inverse-Gamma piror (slightly informative) - the default is non-informative prior
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain,family="gaussian",data=data,verbose=T,mev=vofd,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior1)
#prior2 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000))) # parameter expanded prior (close to non-informative)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain,family="gaussian",data=data,verbose=T,mev=vofd,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior2)
##test runs - with VCV matrix
#prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
#prior4 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,fix=1))) #parameter expanded prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior4)
#summary(m0.0)
#plot(m0.0)

##long runs
prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
model_null_RI2_full_1 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_RI2_full_2 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_RI2_full_3 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_null_RI2_full_1,model_null_RI2_full_2,model_null_RI2_full_3,file="model_null_RI2_full.RData") 

# RI2 - Strain model for the full data set.

#test run
#m0.0<-MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3s)
#summary(m0.0)

#long runs
model_strain_RI2_full_1 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_RI2_full_2 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_RI2_full_3 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
save(model_strain_RI2_full_1,model_strain_RI2_full_2,model_strain_RI2_full_3,file="model_strain_RI2_full.RData") 

# RI2 - full1 model for the full data set (with diet stat and end dates)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full1_RI2_full_1 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_RI2_full_2 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_RI2_full_3 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full1_RI2_full_1,model_full1_RI2_full_2,model_full1_RI2_full_3,file="model_full1_RI2_full.RData") 

# RI2 - full2 model for the full data set (timing of the diet as categorical variable: lactation included or not)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full2_RI2_full_1 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_RI2_full_2 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_RI2_full_3 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full2_RI2_full_1,model_full2_RI2_full_2,model_full2_RI2_full_3,file="model_full2_RI2_full.RData") 

# RI2 - Null model for the chow data subset.

#Define chow data subset variables.
d <- data_NC$RelIntake2_Hd
vofd <- data_NC$RelIntake2_Hd_Var_adj
m_NC <- Bodyweight_matrix_NC
AinvG_NC <- solve(m_NC)
AnivG_NC <- as(AinvG_NC,"dgCMatrix")

##test runs - with VCV matrix
#prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
#prior4 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,fix=1))) #parameter expanded prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior4)
#summary(m0.0)
#plot(m0.0)

##long runs
prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
model_null_RI2_chow_1 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_RI2_chow_2 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_RI2_chow_3 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_null_RI2_chow_1,model_null_RI2_chow_2,model_null_RI2_chow_3,file="model_null_RI2_chow.RData") 

# RI2 - Strain model for the chow data subset.

##test runs - with VCV matrix
#m0.0 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3s)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
#prior4 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,fix=1))) #parameter expanded prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior4)
#summary(m0.0)

##long runs
model_strain_RI2_chow_1 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_RI2_chow_2 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_RI2_chow_3 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
save(model_strain_RI2_chow_1,model_strain_RI2_chow_2,model_strain_RI2_chow_3,file="model_strain_RI2_chow.RData") 

# RI2 - full1 model for the chow data subset (with diet stat and end dates)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full1_RI2_chow_1 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_RI2_chow_2 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_RI2_chow_3 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full1_RI2_chow_1,model_full1_RI2_chow_2,model_full1_RI2_chow_3,file="model_full1_RI2_chow.RData") 

# RI2 - full2 model for the chow data subset (timing of the diet as categorical variable: lactation included or not)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full2_RI2_chow_1 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_RI2_chow_2 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_RI2_chow_3 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_Intake_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full2_RI2_chow_1,model_full2_RI2_chow_2,model_full2_RI2_chow_3,file="model_full2_RI2_chow.RData") 


### Offspring body weight - Bodyweight (BW)
#-----------------------------------------------------
  
#  Define full dataset variables.
d <- data$Bodyweight_Hd
vofd <- data$Bodyweight_Hd_Var_adj
m <- Bodyweight_matrix
AinvG <- solve(m)
AnivG <- as(AinvG,"dgCMatrix")

# BW - Null model for the full data set.

##test runs - without VCV matrix
#m0.0<-MCMCglmm(d~1,random=~Study_ID+Strain,family="gaussian",data=data,verbose=T,mev=vofd,nitt=200000,thin=25,burnin=25000,pr=T) #use stronger prior!
#prior1 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002))) #inverse-Gamma piror (slightly informative) - the default is non-informative prior
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain,family="gaussian",data=data,verbose=T,mev=vofd,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior1)
#prior2 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000))) # parameter expanded prior (close to non-informative)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain,family="gaussian",data=data,verbose=T,mev=vofd,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior2)
##test runs - with VCV matrix
#prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
#prior4 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,fix=1))) #parameter expanded prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior4)
#summary(m0.0)
#plot(m0.0)

##long runs
prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
model_null_BW_full_1 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_BW_full_2 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_BW_full_3 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_null_BW_full_1,model_null_BW_full_2,model_null_BW_full_3,file="model_null_BW_full.RData") 

# BW - Strain model for the full data set.

#test run
#m0.0<-MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3s)
#summary(m0.0)

#long runs
model_strain_BW_full_1 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_BW_full_2 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_BW_full_3 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
save(model_strain_BW_full_1,model_strain_BW_full_2,model_strain_BW_full_3,file="model_strain_BW_full.RData") 

# BW - full1 model for the full data set (with diet stat and end dates)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full1_BW_full_1 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_BW_full_2 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_BW_full_3 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full1_BW_full_1,model_full1_BW_full_2,model_full1_BW_full_3,file="model_full1_BW_full.RData") 

# BW - full2 model for the full data set (timing of the diet as categorical variable: lactation included or not)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full2_BW_full_1 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_BW_full_2 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_BW_full_3 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_choice_diet)+factor(Dam_diet_lactation_incl)+Z_Offspr_diet_E+Z_Offspr_diet_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full2_BW_full_1,model_full2_BW_full_2,model_full2_BW_full_3,file="model_full2_BW_full.RData") 

# BW - Null model for the chow data subset.

#Define chow data subset variables.

d <- data_NC$Bodyweight_Hd
vofd <- data_NC$Bodyweight_Hd_Var_adj
m_NC <- Bodyweight_matrix_NC
AinvG_NC <- solve(m_NC)
AnivG_NC <- as(AinvG_NC,"dgCMatrix")

# BW - null model for the chow data subset.

##test runs - with VCV matrix
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
#prior4 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G2=list(V=1,nu=1,alpha.mu=0,alpha.V=1000),G3=list(V=1,fix=1))) #parameter expanded prior with AnivG
#m0.0 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data,ginverse=list(ES_ID=AnivG), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior4)
#summary(m0.0)
#plot(m0.0)

##long runs
prior3 <- list(R=list(V=1,nu=.002),G=list(G1=list(V=1,nu=.002),G2=list(V=1,nu=.002),G3=list(V=1,fix=1))) #inverse-Gamma prior with AnivG
model_null_BW_chow_1 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_BW_chow_2 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_null_BW_chow_3 <- MCMCglmm(d~1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_null_BW_chow_1,model_null_BW_chow_2,model_null_BW_chow_3,file="model_null_BW_chow.RData") 

# BW - Strain model for the chow data subset.

##test runs - with VCV matrix
#m0.0 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3s)
#summary(m0.0)

##long runs
model_strain_BW_chow_1 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_BW_chow_2 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
model_strain_BW_chow_3 <- MCMCglmm(d~Strain-1,random=~Study_ID+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3s)
save(model_strain_BW_chow_1,model_strain_BW_chow_2,model_strain_BW_chow_3,file="model_strain_BW_chow.RData") 

# BW - full1 model for the chow data subset (with diet stat and end dates)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full1_BW_chow_1 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_BW_chow_2 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full1_BW_chow_3 <- MCMCglmm(d~factor(Offspr_sex)+Z_Dam_diet_start+Z_Dam_diet_end+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full1_BW_chow_1,model_full1_BW_chow_2,model_full1_BW_chow_3,file="model_full1_BW_chow.RData") 

# BW - full2 model for the chow data subset (timing of the diet as categorical variable: lactation included or not)

#test run
#m0.0 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=44000,thin=40,burnin=4000,pr=T,prior=prior3)
#summary(m0.0)

#long runs
model_full2_BW_chow_1 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_BW_chow_2 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
model_full2_BW_chow_3 <- MCMCglmm(d~factor(Offspr_sex)+factor(Dam_diet_lactation_incl)+Z_Dam_diet_exp_E+Z_Dam_diet_exp_PNP_kcal+Z_TBodyweight_age_dPC-1,random=~Study_ID+Strain+ES_ID,family="gaussian",data=data_NC,ginverse=list(ES_ID=AnivG_NC), verbose=T,nitt=1100000,thin=1000,burnin=100000,pr=T,prior=prior3)
save(model_full2_BW_chow_1,model_full2_BW_chow_2,model_full2_BW_chow_3,file="model_full2_BW_chow.RData") 
