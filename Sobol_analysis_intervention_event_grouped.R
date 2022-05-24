# This script performs sobol sensitivity analysis of simple TB models
#   uses the sensobol package to generate the matrices (using sobol sequences)
#   uses the sensitivity package to calculate the indices (using the jansen method)

# Looks at the sensitivity of impact of a hypothetical mass screening to: 
# model choice, natural history parameters, input data, intervention parameters
#
# In this version inputs are grouped (model, inputs, intervention parameters, natural history parameters) - uses approach described in:
# https://click.endnote.com/viewer?doi=10.1029%2F2018wr023403&token=WzI5NzA3MzcsIjEwLjEwMjkvMjAxOHdyMDIzNDAzIl0.QX09aHOi2mW-04RpVOwKldk0Zik

# Set the working directory
setwd("~/GitHub/TB_structure_analysis")

# USER DEFINED INPUTS ##########################################################
# Set the number of samples for the group distributions - this is the number of realizations we create for each group
N_group <- 10000
# Set the number of samples for the sobol analysis
N_sobol <- 10000

# Set the data set to use for baseline incidence. Currently one of: global, PHI
inc_source <- "PHI"
# and the distribution to assume. Currently one of: unif, lnorm
inc_dist <- "unif"

# set the time to run the intervention for 
tend <- 10
times=seq(0,tend)
################################################################################

# load libraries
require(deSolve)
require(ggplot2)
require(reshape2)
require(stats)
library(FME)
library(rriskDistributions)
library(epiR)
library(sensobol)
library(sensitivity)
library(extraDistr)
library(dplyr)

# define palette for plotting
colorBlind   <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load the TB models
source("Model_1_intervention_event.R")
source("Model_2_intervention_event.R")
# load the event function for the intervention
source("intervention_event.R")

# fit distributions to parameters
# this takes 95% CIs for parameters from literature and fits distributions to be used in sampling
source("Par_dist_gen.R")

# define parameters in each group
nat_params <- c("RR of re-infection","TB mortality","Self-cure",                          # Natural history
                "Slow progression","Transition to remote infection",
                "Proportion fast","Fast progression (1)","Fast progression (2)")                                        # natural history parameters
ACF_params <- c("ACF sensitivity","Treatment uptake")                                     # ACF 
PT_params <- c("TST completion","TST sensitivity","PT uptake",                            # PT 
               "PT completion","PT efficacy")   
in_params <- c( "Baseline incidence","Treatment success","Case detection ratio")                                                   # external inputs

# now create realizations for the model choice and for each parameter group 

# for model choice create matrix with all possible combinations (4)
mod_mat <- cbind(c(1,1,2,2),c(1,2,1,2))
colnames(mod_mat) <- c("Progression model","Reinfection model")

# for each parameter group create matrices with N_group realizations - use sobol sequences
# sobol matrices are U[0,1] so convert to correct distributions using distributions from Par_dist_gen.R

nat_mat <- sobol_matrices(matrices="A",N=N_group,params=nat_params)
nat_mat[,"RR of re-infection"] <- qbeta(nat_mat[,"RR of re-infection"],fit_q[1],fit_q[2])
nat_mat[,"TB mortality"] <- qlnorm(nat_mat[,"TB mortality"],fit_m[1],fit_m[2])
nat_mat[,"Self-cure"] <- qlnorm(nat_mat[,"Self-cure"],fit_w[1],fit_w[2])
nat_mat[,"Slow progression"] <- qlnorm(nat_mat[,"Slow progression"],fit_v[1],fit_v[2])
nat_mat[,"Transition to remote infection"] <- qlnorm(nat_mat[,"Transition to remote infection"],fit_k[1],fit_k[2])
nat_mat[,"Proportion fast"] <- qbeta(nat_mat[,"Proportion fast"],fit_g[1],fit_g[2])
nat_mat[,"Fast progression (1)"] <- qlnorm(nat_mat[,"Fast progression (1)"],fit_e[1,1],fit_e[2,1])
nat_mat[,"Fast progression (2)"] <- qlnorm(nat_mat[,"Fast progression (2)"],fit_e[1,2],fit_e[2,2])

ACF_mat <- sobol_matrices(matrices="A",N=N_group,params=ACF_params)
ACF_mat[,"ACF sensitivity"] <- qbeta(ACF_mat[,"ACF sensitivity"],fit_ACF_sens[1],fit_ACF_sens[2]) 
ACF_mat[,"Treatment uptake"] <- qbeta(ACF_mat[,"Treatment uptake"],fit_treat_uptake[1],fit_treat_uptake[2]) 

PT_mat <- sobol_matrices(matrices="A",N=N_group,params=PT_params)
PT_mat[,"TST completion"] <- qbeta(PT_mat[,"TST completion"],fit_TST_comp[1],fit_TST_comp[2])
PT_mat[,"TST sensitivity"] <- qbeta(PT_mat[,"TST sensitivity"],fit_TST_sens[1],fit_TST_sens[2]) 
PT_mat[,"PT uptake"] <- qbeta(PT_mat[,"PT uptake"],fit_PT_start[1],fit_PT_start[2]) 
PT_mat[,"PT completion"] <- qbeta(PT_mat[,"PT completion"],fit_PT_comp[1],fit_PT_comp[2]) 
PT_mat[,"PT efficacy"] <- 1-qbeta(PT_mat[,"PT efficacy"],fit_PT_eff[1],fit_PT_eff[2]) 

in_mat <- sobol_matrices(matrices="A",N=N_group,params=in_params)
assign("fit_inc",get(paste("inc_",inc_source,"_",inc_dist,sep="")))
if (inc_dist=="unif")
  in_mat[,"Baseline incidence"] <- qunif(in_mat[,"Baseline incidence"],fit_inc[1],fit_inc[2])/100000
if (inc_dist=="lnorm")
  in_mat[,"Baseline incidence"] <- qlnorm(in_mat[,"Baseline incidence"],fit_inc[1],fit_inc[2])/100000
in_mat[,"Case detection ratio"] <- qunif(in_mat[,"Case detection ratio"],fit_CDR[1],fit_CDR[2]) 
in_mat[,"Treatment success"] <- qunif(in_mat[,"Treatment success"],fit_tau[1],fit_tau[2])

# Now create Sobol matrices for sampling from the realizations for each group
params <- c("Parameters","ACF","PT","Data","Model")
kpars <- length(params)
k2pars <- kpars+2
mat <- sobol_matrices(N = N_sobol,params = params)
# sobol matrices are U[0,1] so convert to uniform integers between 1 and max number of realizations
mat[,"Parameters"] <- qdunif(mat[,"Parameters"],min=1,max=dim(nat_mat)[1])
mat[,"ACF"] <- qdunif(mat[,"ACF"],min=1,max=dim(ACF_mat)[1])
mat[,"PT"] <- qdunif(mat[,"PT"],min=1,max=dim(PT_mat)[1])
mat[,"Data"] <- qdunif(mat[,"Data"],min=1,max=dim(in_mat)[1])
mat[,"Model"] <- qdunif(mat[,"Model"],min=1,max=dim(mod_mat)[1])

# split out the A and B matrices for the sensitivity package
x1 <- mat[1:N_sobol,]
x2 <- mat[(N_sobol+1):(2*N_sobol),]

xin <- soboljansen(model = NULL, x1, x2, nboot = 1000)

# array to store model outputs
y <- mat.or.vec(dim(mat)[1],6)
beta_out <- rep(0,dim(mat)[1])

# Now create matrix x which takes the correct rows of the parameter matrices based on the group matrix
x <- cbind(nat_mat[mat[,"Parameters"],],
           ACF_mat[mat[,"ACF"],],
           PT_mat[mat[,"PT"],],
           in_mat[mat[,"Data"],],
           mod_mat[mat[,"Model"],])
mm <- x[,"Progression model"]
zz <- x[,"Reinfection model"]

# set the intervention coverage
coverage_ACF <- 0.7
coverage_PT <- 0.7

# Run the model for each row of x, return % reduction in TB incidence and mortality (per 100k) after 1 and 10 yrs
source("Run_model.R")

# calculate sobol indices
xtemp <- xin
tell(xtemp,inc_drop_1,return.var=c("S.boot","T.boot"))
x_sens_inc_1g <- xtemp

xtemp <- xin
tell(xtemp,inc_drop_10,return.var=c("S.boot","T.boot"))
x_sens_inc_10g <- xtemp

xtemp <- xin
tell(xtemp,mort_drop_1,return.var=c("S.boot","T.boot"))
x_sens_mort_1g <- xtemp

xtemp <- xin
tell(xtemp,mort_drop_10,return.var=c("S.boot","T.boot"))
x_sens_mort_10g <- xtemp

# calculate indices for a dummy parameter and bootstrap
# use approach in https://click.endnote.com/viewer?doi=10.1016%2Fj.envsoft.2017.02.001&token=WzI5NzA3MzcsIjEwLjEwMTYvai5lbnZzb2Z0LjIwMTcuMDIuMDAxIl0.qaCejmzKjLqhkWIJ3UhiKwdr7qU
# uses model outputs that correspond to input matrices x1 and x2

N_boot <- 1000
Sd_boot <- mat.or.vec(N_boot,4)
Td_boot <- Sd_boot
#combine the 4 outputs of interest so we can do this in a loop
outputs <- cbind(inc_drop_1,inc_drop_10,mort_drop_1,mort_drop_10)

for(i in 1:N_boot){
  boot_i <- sample(seq(1,N_sobol),N_sobol,replace=TRUE)
  for (o in 1:4){
    #calculate f02 - eqn 7
    f02 <- (1/N_sobol)*sum(outputs[boot_i,o]*outputs[(boot_i+N_sobol),o])
    #calculate total variance - eqn 6
    Vtot <- (1/((2*N_sobol)-1))*sum(outputs[boot_i,o]^2+outputs[(boot_i+N_sobol),o]^2)-f02 
    #calculate partial variance Vd - eqn 12
    Vd <- (1/(N_sobol-1))*sum(outputs[boot_i,o]*outputs[(boot_i+N_sobol),o])-f02
    #calculate partial variance V~d - eqn 13
    Vnotd <- (1/(N_sobol-1))*sum(outputs[(boot_i+N_sobol),o]*outputs[(boot_i+N_sobol),o])-f02
    #calculate Sd - eqn 3
    Sd_boot[i,o] <- Vd/Vtot
    #calculate Td - eqn 4
    Td_boot[i,o] <- 1-(Vnotd/Vtot)
  }
}

Sd <- cbind(apply(Sd_boot, 2, quantile, probs = c(0.5),  na.rm = TRUE),
            apply(Sd_boot, 2, quantile, probs = c(0.025),  na.rm = TRUE),
            apply(Sd_boot, 2, quantile, probs = c(0.975),  na.rm = TRUE))

Td <- cbind(apply(Td_boot, 2, quantile, probs = c(0.5),  na.rm = TRUE),
            apply(Td_boot, 2, quantile, probs = c(0.025),  na.rm = TRUE),
            apply(Td_boot, 2, quantile, probs = c(0.975),  na.rm = TRUE))

# combine these so we can add to the plot
dummy_ind <- cbind(rbind(Sd,Td),
                   rep(c(rep("Incidence",2),rep("Mortality",2)),2),
                   c(rep("Si",4),rep("Ti",4)),
                   rep(c("1yr","10yr"),4))
colnames(dummy_ind) <- c("original","min","max","variable","index","time")
dummy_ind <- as.data.frame(dummy_ind)

# combine results for plotting
df1 <- data.frame(cbind(x_sens_inc_1g$S[,c(1,4,5)],params,"Incidence","1yr","Si",0))
df2 <- data.frame(cbind(x_sens_inc_10g$S[,c(1,4,5)],params,"Incidence","10yr","Si",0))
df3 <- data.frame(cbind(x_sens_mort_1g$S[,c(1,4,5)],params,"Mortality","1yr","Si",0))
df4 <- data.frame(cbind(x_sens_mort_10g$S[,c(1,4,5)],params,"Mortality","10yr","Si",0))
df5 <- data.frame(cbind(x_sens_inc_1g$T[,c(1,4,5)],params,"Incidence","1yr","Ti",0))
df6 <- data.frame(cbind(x_sens_inc_10g$T[,c(1,4,5)],params,"Incidence","10yr","Ti",0))
df7 <- data.frame(cbind(x_sens_mort_1g$T[,c(1,4,5)],params,"Mortality","1yr","Ti",0))
df8 <- data.frame(cbind(x_sens_mort_10g$T[,c(1,4,5)],params,"Mortality","10yr","Ti",0))
colnames(df1) <- c("original","min","max","parameter","variable","time","index","sig")               
names(df8) <- names(df7) <- names(df6) <- names(df5) <- names(df4) <- names(df3) <- names(df2) <- names(df1)
ind_to_plot <- rbind(df1,df2,df3,df4,df5,df6,df7,df8)

# Compare input indices to the dummy values
# If they overlap then conclude that input is not important
# For each input and each output, check if the lower value of the Si is less than the upper value of the dummy Si

for (ii in c("Incidence","Mortality")){
  for (jj in c("1yr","10yr")){
    for (kk in c("Si","Ti")){
      
      temp_ind <- ind_to_plot[ind_to_plot$variable==ii&
                                ind_to_plot$time==jj&
                                ind_to_plot$index==kk,]
      
      temp_dummy <-dummy_ind[dummy_ind$variable==ii&
                               dummy_ind$time==jj&
                               dummy_ind$index==kk,]
      
      ind_to_plot[ind_to_plot$variable==ii&
                    ind_to_plot$time==jj&
                    ind_to_plot$index==kk,"sig"] <- as.numeric(as.character(temp_ind$min)) > as.numeric(as.character(temp_dummy$max))
      
    }
  }
}

ind_to_plot$par_order <- factor(ind_to_plot$parameter,
                                levels = c("Model","Parameters","Data","ACF","PT"))

ind_to_plot$time_order <- factor(ind_to_plot$time,
                                 levels=c("1yr","10yr"))

ind_plot <- ggplot(ind_to_plot,aes(par_order,original))+
  geom_col(aes(fill=index,group=index,alpha=as.factor(sig)),colour="black",position=position_dodge())+
  geom_errorbar(aes(ymin=min, ymax=max,group=index), width=.2,position=position_dodge(0.9))+ 
  facet_grid(time_order~variable)+
  scale_x_discrete(guide = guide_axis(angle = 90))+
  ylab("Value of index")+
  xlab("")+
  scale_fill_manual(values=colorBlind,name="")+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(strip.background =element_rect(fill="white"))+
  guides(alpha = "none")+
  scale_y_continuous(expand = c(0, 0))+
  scale_alpha_manual(values=c(0.2,0.8))

# Check convergence of sobol indices by recalculating for sub-samples of the original matrix
# Due to sampling scheme, need to select appropriate elements for sub-sample (i.e. can;t just pick randomly)
sub_step <- 100
sub_seq <- seq(sub_step,N_sobol,sub_step)  # sequence of sub-sample sizes to use

sub_S_inc_1 <- array(0,c(length(sub_seq),kpars,3))
sub_S_inc_10 <- sub_S_inc_1
sub_S_mort_1 <- sub_S_inc_1
sub_S_mort_10 <- sub_S_inc_1
sub_T_inc_1 <- sub_S_inc_1
sub_T_inc_10 <- sub_S_inc_1
sub_T_mort_1 <- sub_S_inc_1
sub_T_mort_10 <- sub_S_inc_1

for(i in 1:length(sub_seq)){
  
  N_sub <- sub_seq[i]       # set size of sub-sample 
  x1sub <- x1[1:N_sub,]     # get x1
  x2sub <- x2[1:N_sub,]     # get x2
  xsub <- soboljansen(model = NULL, x1sub, x2sub, nboot = 1000) # generate x
  # get the rows of the model output to use for the sub-sample
  tt <- split(seq(1:dim(xin$X)[1]),ceiling(seq_along(seq(1:dim(xin$X)[1]))/N_sobol))
  tt <- matrix(unlist(tt),ncol=k2pars)
  ttt <- as.vector(tt[1:N_sub,])
  # re-calculate indices
  Ysub <- inc_drop_1[ttt]
  xtemp <- xsub
  tell(xtemp,Ysub)
  sub_S_inc_1[i,,1] <- xtemp$S[,"original"]
  sub_T_inc_1[i,,1] <- xtemp$T[,"original"]
  sub_S_inc_1[i,,2] <- xtemp$S[,"min. c.i."]
  sub_T_inc_1[i,,2] <- xtemp$T[,"min. c.i."]
  sub_S_inc_1[i,,3] <- xtemp$S[,"max. c.i."]
  sub_T_inc_1[i,,3] <- xtemp$T[,"max. c.i."]

  Ysub <- inc_drop_10[ttt]
  xtemp <- xsub
  tell(xtemp,Ysub)
  sub_S_inc_10[i,,1] <- xtemp$S[,"original"]
  sub_T_inc_10[i,,1] <- xtemp$T[,"original"]
  sub_S_inc_10[i,,2] <- xtemp$S[,"min. c.i."]
  sub_T_inc_10[i,,2] <- xtemp$T[,"min. c.i."]
  sub_S_inc_10[i,,3] <- xtemp$S[,"max. c.i."]
  sub_T_inc_10[i,,3] <- xtemp$T[,"max. c.i."]
  
  Ysub <- mort_drop_1[ttt]
  xtemp <- xsub
  tell(xtemp,Ysub)
  sub_S_mort_1[i,,1] <- xtemp$S[,"original"]
  sub_T_mort_1[i,,1] <- xtemp$T[,"original"]
  sub_S_mort_1[i,,2] <- xtemp$S[,"min. c.i."]
  sub_T_mort_1[i,,2] <- xtemp$T[,"min. c.i."]
  sub_S_mort_1[i,,3] <- xtemp$S[,"max. c.i."]
  sub_T_mort_1[i,,3] <- xtemp$T[,"max. c.i."]
  
  Ysub <- mort_drop_10[ttt]
  xtemp <- xsub
  tell(xtemp,Ysub)
  sub_S_mort_10[i,,1] <- xtemp$S[,"original"]
  sub_T_mort_10[i,,1] <- xtemp$T[,"original"]
  sub_S_mort_10[i,,2] <- xtemp$S[,"min. c.i."]
  sub_T_mort_10[i,,2] <- xtemp$T[,"min. c.i."]
  sub_S_mort_10[i,,3] <- xtemp$S[,"max. c.i."]
  sub_T_mort_10[i,,3] <- xtemp$T[,"max. c.i."]
  
}

# rearrange sub-sample data to plot
# incidence
tempS <- cbind(sub_seq,sub_S_inc_1[,,1])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS1 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_inc_1[,,2])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS2 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_inc_1[,,3])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS3 <- melt(tempS,id.vars=c("N"))
tempSI1 <- cbind(tempS1,tempS2[,3],tempS3[,3],"Incidence","1yr","Si")
colnames(tempSI1) <- c("N","parameter","original","min","max","variable","time","index")

tempS <- cbind(sub_seq,sub_S_inc_10[,,1])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS1 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_inc_10[,,2])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS2 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_inc_10[,,3])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS3 <- melt(tempS,id.vars=c("N"))
tempSI10 <- cbind(tempS1,tempS2[,3],tempS3[,3],"Incidence","10yr","Si")
colnames(tempSI10) <- c("N","parameter","original","min","max","variable","time","index")

tempT <- cbind(sub_seq,sub_T_inc_1[,,1])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT1 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_inc_1[,,2])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT2 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_inc_1[,,3])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT3 <- melt(tempT,id.vars=c("N"))
tempTI1 <- cbind(tempT1,tempT2[,3],tempT3[,3],"Incidence","1yr","Ti")
colnames(tempTI1) <- c("N","parameter","original","min","max","variable","time","index")

tempT <- cbind(sub_seq,sub_T_inc_10[,,1])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT1 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_inc_10[,,2])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT2 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_inc_10[,,3])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT3 <- melt(tempT,id.vars=c("N"))
tempTI10 <- cbind(tempT1,tempT2[,3],tempT3[,3],"Incidence","10yr","Ti")
colnames(tempTI10) <- c("N","parameter","original","min","max","variable","time","index")

tempI <- rbind(tempSI1,tempSI10,tempTI1,tempTI10)

# Mortality
tempS <- cbind(sub_seq,sub_S_mort_1[,,1])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS1 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_mort_1[,,2])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS2 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_mort_1[,,3])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS3 <- melt(tempS,id.vars=c("N"))
tempSM1 <- cbind(tempS1,tempS2[,3],tempS3[,3],"Mortality","1yr","Si")
colnames(tempSM1) <- c("N","parameter","original","min","max","variable","time","index")

tempS <- cbind(sub_seq,sub_S_mort_10[,,1])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS1 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_mort_10[,,2])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS2 <- melt(tempS,id.vars=c("N"))
tempS <- cbind(sub_seq,sub_S_mort_10[,,3])
colnames(tempS) <- c("N",params)
tempS <- as.data.frame(tempS)
tempS3 <- melt(tempS,id.vars=c("N"))
tempSM10 <- cbind(tempS1,tempS2[,3],tempS3[,3],"Mortality","10yr","Si")
colnames(tempSM10) <- c("N","parameter","original","min","max","variable","time","index")

tempT <- cbind(sub_seq,sub_T_mort_1[,,1])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT1 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_mort_1[,,2])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT2 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_mort_1[,,3])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT3 <- melt(tempT,id.vars=c("N"))
tempTM1 <- cbind(tempT1,tempT2[,3],tempT3[,3],"Mortality","1yr","Ti")
colnames(tempTM1) <- c("N","parameter","original","min","max","variable","time","index")

tempT <- cbind(sub_seq,sub_T_mort_10[,,1])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT1 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_mort_10[,,2])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT2 <- melt(tempT,id.vars=c("N"))
tempT <- cbind(sub_seq,sub_T_mort_10[,,3])
colnames(tempT) <- c("N",params)
tempT <- as.data.frame(tempT)
tempT3 <- melt(tempT,id.vars=c("N"))
tempTM10 <- cbind(tempT1,tempT2[,3],tempT3[,3],"Mortality","10yr","Ti")
colnames(tempTM10) <- c("N","parameter","original","min","max","variable","time","index")

tempM <- rbind(tempSM1,tempSM10,tempTM1,tempTM10)

temp <- rbind(tempI,tempM)

#temp$N%in%seq(1000,20000,1000)&
tthk <- temp[temp$index=="Si",]
ttyk <- ind_to_plot[ind_to_plot$index=="Si",]

# Plot with lines to show final values and CI

tthk$time_order <- factor(tthk$time,
                                 levels=c("1yr","10yr"))

conv_plot <- ggplot(tthk,aes(N/1000,original,color=parameter))+
  #geom_errorbar(aes(ymin=min,ymax=max,color=parameter))+
  #geom_point(alpha=0.5,shape=20)+
  geom_line()+
  facet_grid(time_order~variable,scale="free")+
  ylab("Value of index")+xlab("N (thousands)")+
  #geom_hline(data=ttyk,aes(yintercept=original,color=parameter),linetype="dashed")+
  #geom_hline(data=ttyk,aes(yintercept=min,color=parameter))+
  #geom_hline(data=ttyk,aes(yintercept=max,color=parameter))+
  scale_color_manual(values=colorBlind)+
  theme_bw()+
  theme(legend.position="bottom")+
  theme(strip.background =element_rect(fill="white"))+
  theme(legend.title = element_blank())+
  scale_x_continuous(limits=c(0,10))
  