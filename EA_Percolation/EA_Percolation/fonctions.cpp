//
//  fonctions.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include "fonctions.hpp"

std::vector<double> tabTau1(int Nterme, int d){
    
    std::vector<double> tab(Nterme);
    tab[Nterme-1] = 1.0/ Nterme;
    for(int i = Nterme - 2; i != -1 ; i--){
        double b = borneInfS(i+1, d);
        tab[i] = (1.0 + tab[i+1]*b)/(b+i+1);
    }
    
    return tab;
}
double Tau1(int Nterme, int d){
    std::vector<double> tab = tabTau1(Nterme, d);
    
    return tab.front();
    
}

double borneInfS(int i, int d){
    return ceil(2*(d-1)*pow(i, (d-2)/(double)(d-1)));
}

std::vector< std::vector<double> > tabTau2(int Niterme, int Njterme, int d ){
    std::vector< std::vector<double> > tab(Niterme, std::vector<double>(Njterme) );
    tab[0][0] = 2.0;
    std::cout << tab[0][0] << std::endl; 
    
    return tab;
}
