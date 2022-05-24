# Calculate the steady state for model 1 - see appendix to paper for derivation of this

# Constants used in the solution
Qq <- w+(d*tau)+u+m
I <- incidence/Qq
Zq <- (p*e+u)*(u+m*I)*e*(1-g)
Yq <- (p*e+u)*(e+u)
Wq <- (e+u)*p*e*q*(1-g)

# Parts of quadratic in beta
cq <- Qq*Yq*u*(v+u) - Yq*v*(u*(w+d*tau)) 

bq <- Qq*Yq*(u*q*(1-g)*I + I*(v+u)) - Zq*(v+u) - Yq*v*(g*(u+m*I) + (w+d*tau)*I) - Wq*(w+d*tau)*u*I  

aq <- Qq*Yq*(q*(1-g)*I*I) - Zq*q*(1-g)*I - Wq*((u+m*I)*g*I + (w+d*tau)*I*I)

# Solution to quadratic for beta
betaq <- max(c(0,(-bq+sqrt(bq*bq - 4*aq*cq))/(2*aq),(-bq-sqrt(bq*bq - 4*aq*cq))/(2*aq)))
# And values of other variables
S <- (u+m*I)/(u+betaq*I)
LF <- (1-g)*betaq*I*S/(e+u)
LS <- (g*betaq*I*S + (w+d*tau)*I)/(q*(1-g)*betaq*I + v+ u)
LR <- q*(1-g)*betaq*I*LS/(p*e+u)

I <- incidence/(w+d+u+m)
It <- I*d*(1-tau)/(d*tau+w+u+m)



