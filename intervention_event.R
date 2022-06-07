# define the intervention event function 
event_func <- function(t, y, parms){
  with (as.list(y),{
    
    ACF_I_success <- I*ACF*tau
    ACF_It_success <- It*ACF*tau
    ACF_I_fail <- I*ACF*(1-tau)
    
    PT_LF <- PT*LF
    PT_LS <- PT*LS
    PT_LR <- PT*LR
    
    S <- S
    LF <- LF - PT_LF
    LS <- LS - PT_LS + ACF_I_success + ACF_It_success
    LR <- LR - PT_LR
    I <- I - ACF_I_success - ACF_I_fail
    It <- It + ACF_I_fail - ACF_It_success
    C <- C + PT_LF + PT_LS + PT_LR
    
    return(c(S,LF,LS,LR,I,It,C))
  })
}