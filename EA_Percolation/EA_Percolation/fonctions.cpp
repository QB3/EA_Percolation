//
//  fonctions.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include "fonctions.hpp"

//renvoie un tableau de taille Nterme +1
std::vector<double> tabTau1(int Nterme, int d){
    std::vector<double> tab(Nterme+1);
    tab[Nterme] = 1.0/ (Nterme+1);
    for(int i = Nterme - 1; i != 0 ; i--){
        double b = borneInfS(i, d);
        tab[i] = (1.0 + tab[i+1]*b)/(b+i);
    }
    return tab;
}

double Tau1(int Nterme, int d){
    std::vector<double> tab = tabTau1(Nterme, d);
    return tab[1];
}

double borneInfS(int i, int d){
    return ceil(2*(d-1)*pow(i, (d-2)/(double)(d-1)));
}

std::vector< std::vector<double> > tabTau2(int Niterme, int Njterme, int d ){

    std::vector< std::vector<double> > tab(Niterme, std::vector<double>(Njterme+1));
    tab[Niterme-1] = tabTau1(Niterme, d);
    double A=0.787064;
    tab[Niterme-1][0]=sqrt(3 *M_PI)/sqrt(Njterme)+3/(Njterme*A)/exp(-Njterme*A*A/3);

    for(int i = Niterme - 2; i != -1 ; i--){
	tab[i][Njterme]=1.0/Njterme;
	for(int j = Njterme-2; j!=-1; j--){
		double s_i = borneInfS(i+1, d);
		double s_j = borneInfS(j, d);
		double B = std::max(i-j, 0);
		tab[i][j]=(1+s_i * tab[i+1][j] + (s_j + B)*tab[i][j+1])/(s_i+s_j+B+j);
	}
    }
    
    return tab;
}

double Tau2(int Niterme, int Njterme, int d ){
	std::vector< std::vector<double> > tab(tabTau2(Niterme, Njterme, d ));
	return tab[0][0];
}
