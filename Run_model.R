# Runs the model for each row in an input matrix x
# Returns the intervention impacts

for (ii in 1:dim(x)[1]){
  
  # fix background life expectancy
  u <- 1/50

  # assign the core model parameters
  m <- as.numeric(x[ii,"TB mortality"])
  w <- as.numeric(x[ii,"Self-cure"])
  v <- as.numeric(x[ii,"Slow progression"])
  k <- as.numeric(x[ii,"Transition to remote infection"])
  g <- as.numeric(x[ii,"Proportion fast"])
  tau <- as.numeric(x[ii,"Treatment success"])
  
  # assign the inputs
  incidence <- as.numeric(x[ii,"Baseline incidence"])
  CDR <- as.numeric(x[ii,"Case detection ratio"])
  # use CDR, mortality and self cure rates to get the diagnostic rate 
  d <- CDR*(m+w+u)/(1-CDR)
  
  # assign intervention parameters
  ACF <- coverage_ACF*as.numeric(x[ii,"ACF sensitivity"])*
    as.numeric(x[ii,"Treatment uptake"])
  PT <- coverage_PT*as.numeric(x[ii,"TST completion"])*
    as.numeric(x[ii,"TST sensitivity"])*as.numeric(x[ii,"PT uptake"])*
    as.numeric(x[ii,"PT completion"])*as.numeric(x[ii,"PT efficacy"])
  
  # set parameters for protection depending on the model type
  if (zz[ii]==1){
    q <- as.numeric(x[ii,"RR of re-infection"])
    p <- 1
  }
  if (zz[ii]==2){
    q <- 1
    p <- as.numeric(x[ii,"RR of re-infection"])
  }
  
  # set the parameters for progression and run the model - depending on the model type
  if (mm[ii]==1){
    e <- as.numeric(x[ii,"Fast progression (1)"])
    source("SS_1_beta.R")      # find the steady state by solving for beta (no intervention)
    beta_out[ii] <- betaq
    pars <- c(beta = betaq)    # pass beta to the parameters               
    state <- c(S=S,LF=LF,LS=LS,LR=LR,I=I,It=It,C=0) # set the ICs based on the steady state
    
    temp <- ode(func = derivs_1, y = state, times = times, parms = pars,
                events = list(func = event_func, time = 0)) # run the model with the intervention
    
  } 
  if (mm[ii]==2) {
    e <- as.numeric(x[ii,"Fast progression (2)"])
    source("SS_2_beta.R")      # find the steady state by solving for beta (no intervention)
    beta_out[ii] <- betaq
    pars <- c(beta = betaq)    # pass beta to the parameters
    state <- c(S=S,LF=LF,LS=LS,LR=LR,I=I,It=It,C=0) # set the ICs based on the steady state
    
    temp <- ode(func = derivs_2, y = state, times = times, parms = pars,
                events = list(func = event_func, time = 0)) # run the model with the intervention
  }
  
  # Store the output we want
  y[ii,] <- c(temp[1,"Inc"],temp[2,"Inc"],temp[(tend+1),"Inc"],temp[1,"Mort"],temp[2,"Mort"],temp[(tend+1),"Mort"]) 
}

# name outputs
colnames(y) <- c("Inc_0","Inc_1","Inc_10","Mort_0","Mort_1","Mort_10")

# calculate % change in incidence over 1 year and 10 years
inc_drop_1 <- 100*(y[,"Inc_0"]-y[,"Inc_1"])/y[,"Inc_0"]
inc_drop_10 <- 100*(y[,"Inc_0"]-y[,"Inc_10"])/y[,"Inc_0"]
# calculate % change in mortality over 1 year and 10 years
mort_drop_1 <- 100*(y[,"Mort_0"]-y[,"Mort_1"])/y[,"Mort_0"]
mort_drop_10 <- 100*(y[,"Mort_0"]-y[,"Mort_10"])/y[,"Mort_0"]