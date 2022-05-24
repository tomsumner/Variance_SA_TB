# Model 2 

derivs_2 <- function(t,state,pars){
  with(as.list(c(state,pars)),{
  
    # Check total population stays at 1
    Total = S+LF+LS+LR+I+It+C
      
    # Births set to total deaths to keep total population constant 
    births <- u*Total + m*(I+It) 
      
    # force of infection
    FOI <- beta*(I+It)
      
    # Susceptible
    dS <- births - u*S -                    # plus births minus deaths 
            FOI*S                             # minus infection
      
    # Latent early
    dLF <- (1-g)*S*FOI -                    # plus infection from S
           e*LF -                           # minus progression from LA to disease
           u*LF                             # minus death
      
    # Latent late
    dLS <- g*S*FOI +                        # plus infection from S
           g*q*C*FOI -                      # plus reinfection from C
           v*LS -                           # minus progression from LB to disease    
           q*(1-g)*LS*FOI -                 # minus reinfection 
           u*LS +                           # minus death  
           w*(I+It) + d*tau*(I+It)          # plus self-cure from I and successful treatment 
             
    # Latent reinfection
    dLR <- q*(1-g)*LS*FOI +                 # plus reinfection from LS
           q*(1-g)*C*FOI -                  # plus reinfection from C
           p*e*LR -                         # minus progression to disease
           u*LR                             # minus death
      
    # Infectious
    dI <- e*LF + v*LS + p*e*LR -            # plus progression from LA and LB
          w*I - d*I - u*I - m*I             # minus self-cure, diagnosis and death
        
    # Infectious failed treatment 
    dIt <- d*(1-tau)*I -                    # plus failed treatment from I
           w*It - d*tau*It - u*It - m*It    # minus self-cure, diagnosis (and effective treatment) and death  
      
    # Post preventive therapy
    dC <- -q*C*FOI -                        # minus reinfection
          u*C                               # minus death
      
    # Derived outputs
    Inc <- e*LF + v*LS + p*e*LR             # Total Incidence
    Notif <- d*I                            # Notifications
    Mort <- m*(I+It)                        # Mortality
      
    # Return these things
    return(list(c(dS, dLF, dLS, dLR, dI, dIt, dC), 
                Total = Total,
                Inc = Inc,
                Notif = Notif,
                Mort = Mort))
  })
}
  
