# Calculate the steady state for model 1 - see appendix to paper for derivation of this

# Constants used in the solution
Zq <- k+(p*e)+u
Yq <- k+e+u
Qq <- w+(d*tau)+u+m
I <- incidence/Qq

# Parts of quadratic in beta
cq <- Qq*Yq*(Zq*u*(u+v)) - Yq*(w+d*tau)*(Zq*v*u)

bq <- Qq*Yq*(Zq*u*q*I + Zq*I*(v+u)-u*k*I*q) - e*(Zq*u*(v+u)+Zq*m*I*(v+u)) - k*Zq*v*(u+m*I) - Yq*(w+d*tau)*(Zq*v*I+u*p*e*q*I) 

aq <- Qq*Yq*(Zq*q*I*I-k*q*I*I) - e*(Zq*u*q*I + Zq*q*m*I*I - u*k*I*q - m*k*I*I*q) - k*p*e*q*I*(u+m*I) - Yq*(w+d*tau)*p*e*q*I*I

# Solution to quadratic for beta
betaq <- max(c(0,(-bq+sqrt(bq*bq - 4*aq*cq))/(2*aq),(-bq-sqrt(bq*bq - 4*aq*cq))/(2*aq)))
# And values of other variables
S <- (u+m*I)/(u+betaq*I)
LF <- betaq*I*S/(Yq)
LS <- (k*betaq*I*Zq*S + (w+(d*tau))*Yq*Zq*I)/(Yq*Zq*(q*betaq*I+v+u)-Yq*k*betaq*I*q)
LR <- q*betaq*I*LS/Zq

I <- incidence/(w+d+u+m)
It <- I*d*(1-tau)/(d*tau+w+u+m)





