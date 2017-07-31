source("~/EA_Percolation/fonctions.R")


terme=10**3
dimension=c(2,3,4,5,22,25,30,35)
#dimension=c(20)

res=NULL
tau_1=NULL
system.time(
#on regarde la valeur de la borne pour plusieurs valeurs de d différentes
for(d in dimension){
  borne=borne_tau_1(terme, d)
  tau_1=c(tau_1, borne)
  res=c(res, 2*borne*d/log(2*d))
}
)
round(t(rbind(dimension,  0.33133313/sqrt(dimension), tau_1)), 4)

#T=tab_borne_tau_1(terme, 5)
