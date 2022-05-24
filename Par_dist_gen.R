# This script fits distributions to the parameter ranges

library(rriskDistributions)

# CORE MODEL PARAMETERS
##########################################
# RR of re-infection/progression (VYNYYCKY AND FINE) - proportion 0-1, beta
q_data <-  c(0.55511,0.59973,0.66696)
fit_q <- get.beta.par(p=c(0.025,0.5,0.975),q=q_data,show.output=FALSE,plot=FALSE)
q_range <- ((100*(q_data[2]-q_data[1])/q_data[2])+(100*(q_data[3]-q_data[2])/q_data[2]))/2
# TB mortality - rate >0, log-normal
m_data <- c(0.112,0.231,0.35)
fit_m <- get.lnorm.par(p=c(0.025,0.5,0.975),q=m_data,show.output=FALSE,plot=FALSE)
m_range <- ((100*(m_data[2]-m_data[1])/m_data[2])+(100*(m_data[3]-m_data[2])/m_data[2]))/2
# self-cure - rate >0, log-normal
w_data <- c(0.208,0.329,0.45)
fit_w <- get.lnorm.par(p=c(0.025,0.5,0.975),q=w_data,show.output=FALSE,plot=FALSE)
w_range <- ((100*(w_data[2]-w_data[1])/w_data[2])+(100*(w_data[3]-w_data[2])/w_data[2]))/2
# e - fast progression, differs between models - rate >0, log-normal
# # Based on Vynnycky and Fine
# e_data <- cbind(c(0.02604,0.02754,0.02873),
#                 c(0.19,0.2,0.21))
# fit_e <- cbind(get.lnorm.par(p=c(0.025,0.5,0.975),q=e_data[,1],show.output=FALSE,plot=FALSE),
#                get.lnorm.par(p=c(0.025,0.5,0.975),q=e_data[,2],show.output=FALSE,plot=FALSE))
# Based on ragonnet et al 
e_data <- cbind(c(0.3066,0.4015,0.5475),
                c(3.358,4.015,5.475))
fit_e <- cbind(get.lnorm.par(p=c(0.025,0.5,0.975),q=e_data[,1],show.output=FALSE,plot=FALSE),
               get.lnorm.par(p=c(0.025,0.5,0.975),q=e_data[,2],show.output=FALSE,plot=FALSE))
e1_range <- ((100*(e_data[2,1]-e_data[1,1])/e_data[2,1])+(100*(e_data[3,1]-e_data[2,1])/e_data[2,1]))/2
e2_range <- ((100*(e_data[2,2]-e_data[1,2])/e_data[2,2])+(100*(e_data[3,2]-e_data[2,2])/e_data[2,2]))/2
# v - slow progression - rate >0, log-normal
# # Based on Vynnycky and Fine
#v_data <- c(0.00029,0.00029904,0.00031)
#fit_v <- get.lnorm.par(p=c(0.025,0.5,0.975),q=v_data,show.output=FALSE,plot=FALSE)
# Based on ragonnet et al 
v_data <- c(0.000913,0.002008,0.004015)
fit_v <- get.lnorm.par(p=c(0.025,0.5,0.975),q=v_data,show.output=FALSE,plot=FALSE)
v_range <- ((100*(v_data[2]-v_data[1])/v_data[2])+(100*(v_data[3]-v_data[2])/v_data[2]))/2
# # Based on Vynnycky and Fine
# k - transition from fast to slow (model 1 only) - rate >0, log-normal
# k_data <- c(0.17127,0.17246,0.17396)
# fit_k <- get.lnorm.par(p=c(0.025,0.5,0.975),q=k_data,show.output=FALSE,plot=FALSE)
# Based on ragonnet et al 
k_data <- c(3.1025,4.015,5.11)
fit_k <- get.lnorm.par(p=c(0.025,0.5,0.975),q=k_data,show.output=FALSE,plot=FALSE)
k_range <- ((100*(k_data[2]-k_data[1])/k_data[2])+(100*(k_data[3]-k_data[2])/k_data[2]))/2
# g - proportion entering fast (model 2 only) - proportion 0-1, beta
# # Based on Vynnycky and Fine
# g_data <- c(0.85636,0.8623,0.8698)
# fit_g <- get.beta.par(p=c(0.025,0.5,0.975),q=g_data,show.output=FALSE,plot=FALSE)
# Based on ragonnet et al 
g_data <- c(0.89,0.91,0.93)
fit_g <- get.beta.par(p=c(0.025,0.5,0.975),q=g_data,show.output=FALSE,plot=FALSE)
g_range <- ((100*(g_data[2]-g_data[1])/g_data[2])+(100*(g_data[3]-g_data[2])/g_data[2]))/2

# ACF PARAMETERS #############################
# ACF_sensitivity - assuming use Xpert for everyone - from https://apps.who.int/iris/bitstream/handle/10665/340243/9789240022713-eng.pdf
# proportion 0-1, beta
ACF_sens_data <- c(0.48,0.69,0.86)
fit_ACF_sens <- get.beta.par(p=c(0.025,0.5,0.975),q=ACF_sens_data,show.output=FALSE,plot=FALSE)
ACF_sens_range <- ((100*(ACF_sens_data[2]-ACF_sens_data[1])/ACF_sens_data[2])+(100*(ACF_sens_data[3]-ACF_sens_data[2])/ACF_sens_data[2]))/2
# treatment uptake following ACF diagnosis 
# from https://apps.who.int/iris/bitstream/handle/10665/340243/9789240022713-eng.pdf - use Gopi study values (pg 13 of WHO doc)
# 1-loss to follow up values from doc
# proportion 0-1, beta
Treat_uptake_data <-c(0.71,0.77,0.82)
fit_treat_uptake <- get.beta.par(p=c(0.025,0.5,0.975),q=Treat_uptake_data,show.output=FALSE,plot=FALSE)
Treat_uptake_range <- ((100*(Treat_uptake_data[2]-Treat_uptake_data[1])/Treat_uptake_data[2])+(100*(Treat_uptake_data[3]-Treat_uptake_data[2])/Treat_uptake_data[2]))/2
# PT PARAMETERS #############################
# TST_sensitivity - from https://click.endnote.com/viewer?doi=10.7326%2F0003-4819-149-3-200808050-00241&token=WzI5NzA3MzcsIjEwLjczMjYvMDAwMy00ODE5LTE0OS0zLTIwMDgwODA1MC0wMDI0MSJd.MkPoqzrvSrOUrNJ3GeNXSE-_lWo
# proportion 0-1, beta
TST_sens_data <- c(0.71,0.77,0.82)
fit_TST_sens <- get.beta.par(p=c(0.025,0.5,0.975),q=TST_sens_data,show.output=FALSE,plot=FALSE)
TST_sens_range <- ((100*(TST_sens_data[2]-TST_sens_data[1])/TST_sens_data[2])+(100*(TST_sens_data[3]-TST_sens_data[2])/TST_sens_data[2]))/2
# The following are from:
# https://www-sciencedirect-com.ez.lshtm.ac.uk/science/article/pii/S147330991630216X?via%3Dihub
# Values for general population cohorts (see supplementary tables)
# all proportions 0-1, beta
# TST completion (of those eligible) [table 1]
TST_comp_data <- c(0.37,0.621,0.87)
fit_TST_comp <- get.beta.par(p=c(0.025,0.5,0.975),q=TST_comp_data,show.output=FALSE,plot=FALSE)
TST_comp_range <- ((100*(TST_comp_data[2]-TST_comp_data[1])/TST_comp_data[2])+(100*(TST_comp_data[3]-TST_comp_data[2])/TST_comp_data[2]))/2
# PT start (of those recommended) [table S8]
PT_start_data <- c(0.75,0.861,0.97)
fit_PT_start <- get.beta.par(p=c(0.025,0.5,0.975),q=PT_start_data,show.output=FALSE,plot=FALSE)
PT_start_range <- ((100*(PT_start_data[2]-PT_start_data[1])/PT_start_data[2])+(100*(PT_start_data[3]-PT_start_data[2])/PT_start_data[2]))/2
# PT complete (of those started) [table S9]
PT_comp_data <- c(0.38,0.546,0.71)
fit_PT_comp <- get.beta.par(p=c(0.025,0.5,0.975),q=PT_comp_data,show.output=FALSE,plot=FALSE)
PT_comp_range <- ((100*(PT_comp_data[2]-PT_comp_data[1])/PT_comp_data[2])+(100*(PT_comp_data[3]-PT_comp_data[2])/PT_comp_data[2]))/2
# 3HP treatment efficacy from:
#https://click.endnote.com/viewer?doi=10.7326%2Fm17-0609&token=WzI5NzA3MzcsIjEwLjczMjYvbTE3LTA2MDkiXQ.P9luno2qeTZhsP1GCucft8ZWdIk
# This is the OR for TB (vs no treatment for RPT-INH)
# In the model we use 1-OR as the proportion who are effectively treated 
# proportion 0-1, beta
PT_eff_data <- c(0.18,0.36,0.73)
fit_PT_eff <- get.beta.par(p=c(0.025,0.5,0.975),q=PT_eff_data,show.output=FALSE,plot=FALSE)
PT_eff_range <- ((100*(PT_eff_data[2]-PT_eff_data[1])/PT_eff_data[2])+(100*(PT_eff_data[3]-PT_eff_data[2])/PT_eff_data[2]))/2
# INPUT DATA ####################################################
# Incidence /100,000 and assuming uniform range between high and low estimates
inc_global_data <- c(114,127,140)  # global range
inc_global_unif <- inc_global_data[c(1,3)]
inc_global_lnorm <- get.lnorm.par(p=c(0.025,0.5,0.975),q=inc_global_data,show.output=FALSE,plot=FALSE)
inc_global_range <- ((100*(inc_global_data[2]-inc_global_data[1])/inc_global_data[2])+(100*(inc_global_data[3]-inc_global_data[2])/inc_global_data[2]))/2

inc_PHI_data <- c(306,539,838)     # Philipines - example of high incidence setting
inc_PHI_unif <- inc_PHI_data[c(1,3)]
inc_PHI_lnorm <- get.lnorm.par(p=c(0.025,0.5,0.975),q=inc_PHI_data,show.output=FALSE,plot=FALSE)
inc_PHI_range <- ((100*(inc_PHI_data[2]-inc_PHI_data[1])/inc_PHI_data[2])+(100*(inc_PHI_data[3]-inc_PHI_data[2])/inc_PHI_data[2]))/2
# CDR calculated from number of new and relapse notified cases and estimated incident cases
fit_CDR <- c(0.53,0.66) # 5824745 new and relapse notified cases, 9870 (8880-10900) thousand incident cases
fit_CDR_mid <- (fit_CDR[1]+fit_CDR[2])/2 
CDR_range <- 100*((fit_CDR_mid-fit_CDR[1])/fit_CDR_mid)
# tau set to be uniform between 80% and 90%
fit_tau <- c(0.8,0.9)
fit_tau_mid <- (fit_tau[1]+fit_tau[2])/2 
tau_range <- 100*((fit_tau_mid-fit_tau[1])/fit_tau_mid)

inputs_range <- c(0,0,q_range,m_range,w_range,v_range,k_range,g_range,
                  e1_range,e2_range,inc_global_range,tau_range,CDR_range,
                  ACF_sens_range,Treat_uptake_range,TST_comp_range,TST_sens_range,
                  PT_start_range,PT_comp_range,PT_eff_range)
