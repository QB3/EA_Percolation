source("/home/qbe/Bureau/EA_Percolation/fonctions.R")

terme=1000
dimension=c(2,3,4,5,7,10,20,40)

res=NULL
tau_1=NULL

#on regarde la valeur de la borne pour plusieurs valeurs de d différentes
for(d in dimension){
  borne=borne_tau_1(terme, d)
  tau_1=c(tau_1, borne)
  res=c(res, 2*borne*d/log(2*d))
}

t(rbind(dimension, tau_1, res))

T=tab_borne_tau_1(terme, 5)
