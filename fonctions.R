#fonction qui calcule les bornes sur les s_i
borne_inf_S=function(nombre_noeuds_infectes, dimension){
  dim=dimension-1
  i=nombre_noeuds_infectes
  return(ceiling(2*dim*i**((dim-1)/dim)))
}

borne_T_k_0=function(k){
  f=function(x){
    res=exp(-k*x)*(x+1)^k
    return(res)
  }
  integrate(f, 0, Inf)$value
}

#on crée la fonction qui nous calcule la borne sur tau_1
borne_tau_1=function(nb_terme, dim){
  last_terme=1/(nb_terme+1)
  borne=last_terme
  for(i in nb_terme:1){
    s=borne_inf_S(i, dim)
    borne=(1+s*borne)/(s+i)
  }
  return(borne)
}

#fonctions qui renvoie un tableau de taille nb_terme+1
#tab[i] contient un majorant de T_i
tab_borne_tau_1=function(nb_terme, dim){
  tab=matrix(0, 1, nb_terme+1)
  tab[nb_terme+1]=1/(nb_terme+1)
  
  for(i in nb_terme:1){
    s=borne_inf_S(i, dim)
    tab[i]=(1+s*tab[i+1])/(s+i)
  }
  return(tab)
}
# tab_borne_tau_1(100, 5)

#attention borne_T_k_0 au delà de 100, valeur interpolée
borne_tau_2=function(d, n_termes_i, n_termes_j){
  #on calcule une borne fine sur T_k_0
  # borne_fine_T_k_0=borne_T_k_0(n_termes_i +1)
  borne_fine_T_k_0=2/sqrt(n_termes_i+1)
  #on initialise T_n_i_k
  tab_bas=tab_borne_tau_1(n_termes_j-1, d)
  #on prend le minimum
  tab_bas=pmin(tab_bas, borne_fine_T_k_0)
  tab_bas=c(borne_fine_T_k_0, tab_bas)
  tab_haut=tab_bas
  for (i in n_termes_i:1){
    tab_haut[n_termes_j+1]=1/n_termes_j
    for( j in n_termes_j:1){
      s_i=borne_inf_S(i, d)
      s_j=borne_inf_S(j-1, d)
      B=max(i-j+1, 0)
      tab_haut[j]=(1+tab_bas[j]*s_i+(s_j+B)*tab_haut[j+1])/(s_i+s_j+B+j-1)
    }
    tab_bas=tab_haut
  }
  return(tab_bas[1])
}

#renvoie un tableau de taille n_termes_i * (n_termes_j+1)
tab_borne_tau_2=function(d, n_termes_i, n_termes_j){
  borne_fine_T_k_0=2/sqrt(n_termes_i+1)
  
  tab_bas=tab_borne_tau_1(n_termes_j-1, d)
  tab_bas=pmin(tab_bas, borne_fine_T_k_0)
  
  #tab_bas=c(tab_bas[1]+1/(n_termes_i+1), tab_bas)
  tab_bas=c(borne_fine_T_k_0, tab_bas)
  tab_haut=tab_bas
  res=matrix(0, n_termes_i, n_termes_j+1)
  for (i in n_termes_i:1){
    tab_haut[n_termes_j+1]=1/n_termes_j
    for( j in n_termes_j:1){
      s_i=borne_inf_S(i, d)
      s_j=borne_inf_S(j-1, d)
      B=max(i-j+1, 0)
      tab_haut[j]=(1+tab_bas[j]*s_i+(s_j+B)*tab_haut[j+1])/(s_i+s_j+B+j-1)
    }
    res[i,]=tab_haut
    tab_bas=tab_haut
  }
  return(res)
}

Borne_T_l_0_0=function(n, M, alpha){
  x=NULL
  for (i in 1:M){
    x[i]=min(rgamma(n, alpha, 1))
  }
  res=mean(x)
  return(res)
}

#attention n_terme i > 50 impossibles
borne_tau_3=function(d, n_termes_i, n_termes_j, n_termes_k){

   # d=10
   # n_termes_i=1000
   # n_termes_j=1000
   # n_termes_k=1000

  # f=function(x){
  #   res=exp(-(n_termes_i+1)*x)*(x+1+x^2/2)^(n_termes_i+1)
  #   return(res)
  # }
  # borne_T_l_0_0=integrate(f, 0, Inf, stop.on.error=FALSE)$value
  
  borne_T_l_0_0=Borne_T_l_0_0(n_termes_i+1, 10000, 3)
  #on rempli la plan loin d'équation i=n_termes_i
  tab_ref=tab_borne_tau_2(d, n_termes_j, n_termes_k)
  haut_tab=matrix(0,1,n_termes_k)
  haut_tab=c(borne_T_l_0_0, haut_tab)
  tab_ref=rbind(haut_tab, tab_ref)
  
  tab_pres=tab_ref
  tab_loin=pmin(tab_ref, borne_T_l_0_0)
  
  for (i in n_termes_i:1){
    for (j in n_termes_j:1){
      for(k in n_termes_k:1){
        
        s_i=borne_inf_S(i, d)
        s_j=borne_inf_S(j-1, d)
        s_k=borne_inf_S(k-1, d)
        B_0_1=max(i-j+1, 0)
        B_1_2=max(j-k, 0)
        
        tab_pres[j, k]= 1 + s_i * tab_loin[j, k] + (s_j+B_0_1) * tab_pres[j+1, k] + (s_k + B_1_2) * tab_pres[j, k+1]
        tab_pres[j,k]=tab_pres[j, k]/(s_i+s_j+s_k+B_0_1+B_1_2+k-1)
      }
    }
    tab_loin=tab_pres
    print(paste(i, "/", n_termes_i))
  }
  return(tab_pres[1,1])
}


