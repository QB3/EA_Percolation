source("~/EA_Percolation/fonctions.R")

dimension=c(29:30)
res=NULL
tau_3=NULL
n_termes_i=10
n_termes_j=10
n_termes_k=10

#temps pour 100*100*100 : 2min 30 sec


for(d in dimension){
  borne=borne_tau_3(d, n_termes_i, n_termes_j, n_termes_k)
  tau_3=c(tau_3, borne)
  res=c(res, borne*d/log(2*d)*2/3)
}

tau_3=tau_3/3
t(rbind(dimension, tau_3, res))

# 
# ptm <- proc.time()
#   borne=borne_tau_3(30, 1000, 1000, 100)
# proc.time()-ptm
# borne/3
#   
# 0.borne_tau_2(30, 1000, 1000)

# d=10
# n_termes_j=100
# n_termes_i=100
# borne=tab_borne_tau_2(d, n_termes_i, n_termes_j)
# borne[1,1]*d/log(2*d)
# borne[1,1]
