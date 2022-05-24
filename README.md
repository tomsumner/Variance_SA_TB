# Variance_SA_TB

This repository contains the code used in 

"Variance-based sensitivity analysis of TB transmission models", Sumner, T and White, R.G.

Model_1_intervention_event.R <- function for model 1 (serial infected states) with intervention coded as an event (using package FME)

Model_2_intervention_event.R <- fucntion for model 2 (parallel infected states) with intervention coded as an event (using package FME)

SS_1_beta.R <- code to solve steady state solution for model 1 (serial infected states) for beta

SS_2_beta.R <- code to solve steady state solution for model 2 (parallel infected states) for beta

Par_dist_gen <- fits distributions to input quantiles taken from literature, uses package RRiskDistributions 

Run_model <- code to run the model for each row in an input matrix and return the reductions in incidence and mortality over 1 and 10 years

Sobol_analysis_intervention_event_individual.R <- code to calculate the Sobol indicies for individual inputs and generate outputs, uses Sensobol package

Sobol_analysis_intervention_event_grouped.R <- code to calculate the Sobol indicies for grouped inputs and generate outputs, uses Sensobol package


 
