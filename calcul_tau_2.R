source("~/EA_Percolation/fonctions.R")

dimension=c(10, 20, 40)

res=NULL
tau_2_sur_2=NULL
n_termes_j=1000
n_termes_i=1000

for(d in dimension){
  borne=borne_tau_2(d, n_termes_i, n_termes_j)
  tau_2_sur_2=c(tau_2_sur_2, borne)
  res=c(res, borne*d/log(2*d))
}

tau_2_sur_2=tau_2_sur_2/2
borne_diag=0.331333333/sqrt(dimension)

t(rbind(dimension, borne_diag,tau_2_sur_2, res))

# borne_tau_2(30, n_termes_i, n_termes_j)

tab=borne_tau_2(20, 10000, 1000)
tab=as.data.frame(tab)
