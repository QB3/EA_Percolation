//
//  fonctions.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include "fonctions.hpp"

using namespace std;

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

//rnevoie i=un tableau de taille Niterme+1 par Njterme +1
vector< vector<double> > tabTau2(int Niterme, int Njterme, int d ){

    vector< vector<double> > tab(Niterme+1, vector<double>(Njterme+1));
    tab[Niterme] = tabTau1(Niterme, d);
    double A=0.787064;
    double borne_fine=sqrt(3 *M_PI)/sqrt(Niterme)+3.0/(Niterme*A)*exp(-Niterme*A*A/3);
    tab[Niterme][0]=borne_fine;

    //rafinement de l'initialisation
    for(int j = 1; j != Njterme ; j++){
	tab[Niterme][j]=std::min(borne_fine, tab[Niterme][j]);
	}

    for(int i = Niterme - 1; i != 0 ; i--){
	tab[i][Njterme]=1.0/Njterme;

	for(int j = Njterme-1; j!=-1; j--){
		double s_i = borneInfS(i, d);
		double s_j = borneInfS(j, d);
		double B = std::max(i-j, 0);
		tab[i][j]=(1+s_i * tab[i+1][j] + (s_j + B)*tab[i][j+1])/(s_i+s_j+B+j);
	}
    }
    
    return tab;
}

double Tau2(int Niterme, int Njterme, int d ){
	vector< vector<double> > tab(tabTau2(Niterme, Njterme, d ));
	return tab[1][0];
}

//renvoie un tableau de taille Niterme * Njterme * Nkterme
vector< vector< vector<double> > > tabTau3(int Niterme, int Njterme, int Nkterme, int d ){
	vector < vector< vector<double> > > tab(Niterme+1, vector< vector<double> >(Njterme+1, vector<double>(Nkterme+1)));
	tab[Niterme] = tabTau2(Njterme,Nkterme, d);
	//à modifier
	double A=1.00738;
	double borne_fine=pow(12, 1.0/3)*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/12)/(Niterme *A*A/12);

	tab[Niterme][0][0]=borne_fine;
	for(int j = 1; j != Njterme ; j++){
		for (int k=0; k!=Nkterme; k++){
			tab[Niterme][j][k] = min(borne_fine, tab[Niterme][j][k]);
		}
	}
	
	for(int i = Niterme - 1; i != 0 ; i--){

		for(int k =0; k!=Nkterme+1; k++){
			tab[i][Njterme][k]=tab[Niterme][Njterme][k];
		}

		for(int j =1;j!=Njterme+1; j++){
			tab[i][j][Nkterme]=tab[Niterme][j][Nkterme];
		}

		for(int j = Njterme-1; j!=-1; j--){
			for(int k =Nkterme-1; k!=-1; k--){ 
				double s_i = borneInfS(i, d);
				double s_j = borneInfS(j, d);
				double s_k = borneInfS(k, d);
				double B_0_1 = max(i-j, 0);
				double B_1_2 = max(j-k, 0);
				tab[i][j][k]=(1 + s_i * tab[i+1][j][k] + (s_j+B_0_1) * tab[i][j+1][k] + (s_k+B_1_2) * tab[i][j][k+1])/(s_i+s_j+s_k+B_0_1+B_1_2+k);
			}
		}
    	}
	return tab;
}

double Tau3(int Niterme, int Njterme, int Nkterme, int d ){
	vector < vector< vector<double> > > tab(tabTau3(Niterme, Njterme, Nkterme, d ));
	return tab[1][0][0];
}

vector< vector< vector<double> > > tabTau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d ){
	vector < vector< vector<double> > > tab(Niterme+1, vector< vector<double> >(Njterme+1, vector<double>(Nkterme+1)));
	vector< vector<double> > tabOpt(tabTau2(NiTau2, NjTau2, d));
	for( int j = 1; j!= Njterme+1; j++){
		for(int k = 0; k!=Nkterme; k++){
			tab[Niterme][j][k]=tabOpt[j][k];
		}
	}

	double A=1.00738;
	double borne_fine=pow(12, 1.0/3)*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/12)/(Niterme *A*A/12);

	tab[Niterme][0][0]=borne_fine;
	for(int j = 1; j != Njterme ; j++){
		for (int k=0; k!=Nkterme; k++){
			tab[Niterme][j][k] = min(borne_fine, tab[Niterme][j][k]);
		}
	}
	
	for(int i = Niterme - 1; i != 0 ; i--){

		for(int k =0; k!=Nkterme+1; k++){
			tab[i][Njterme][k]=tab[Niterme][Njterme][k];
		}

		for(int j =1;j!=Njterme+1; j++){
			tab[i][j][Nkterme]=tab[Niterme][j][Nkterme];
		}

		for(int j = Njterme-1; j!=-1; j--){
			for(int k =Nkterme-1; k!=-1; k--){ 
				double s_i = borneInfS(i, d);
				double s_j = borneInfS(j, d);
				double s_k = borneInfS(k, d);
				double B_0_1 = max(i-j, 0);
				double B_1_2 = max(j-k, 0);
				tab[i][j][k]=(1 + s_i * tab[i+1][j][k] + (s_j+B_0_1) * tab[i][j+1][k] + (s_k+B_1_2) * tab[i][j][k+1])/(s_i+s_j+s_k+B_0_1+B_1_2+k);
			}
		}
    	}
	return tab;
}

double Tau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d ){
	vector < vector< vector<double> > > tab(tabTau3Opt(NiTau2, NjTau2, Niterme, Njterme, Nkterme, d ));
	return tab[1][0][0];
}
