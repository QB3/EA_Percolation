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
		double B = max(i-j, 0);
		tab[i][j]=(1+s_i * tab[i+1][j] + (s_j + B)*tab[i][j+1])/(s_i+s_j+B+j);
	}
    }
    
    return tab;
}

double Tau2(int Niterme, int Njterme, int d ){
	vector< vector<double> > tab(tabTau2(Niterme, Njterme, d ));
	return tab[1][0];
}

//renvoie un tableau de taille Niterme +1 * Njterme +1 * Nkterme + 1
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

//renvoie un tableau de taille Niterme + 1 * Njterme + 1 * Nkterme +1 * Nlterme + 1 
vector< vector< vector<double> > > tabTau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d ){
	vector< vector<double> > tabOpt(tabTau2(NiTau2, NjTau2, d));
	vector < vector< vector<double> > > tab(Niterme+1, vector< vector<double> >(Njterme+1, vector<double>(Nkterme+1)));
	for( int j = 1; j!= Njterme+1; j++){
		for(int k = 0; k!=Nkterme+1; k++){
			tab[Niterme][j][k]=tabOpt[j][k];
		}
	}

	double A=1.00738;
	double borne_fine=pow(12, 1.0/3)*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/12)/(Niterme *A*A/12);

	tab[Niterme][0][0]=borne_fine;
	for(int j = 1; j != Njterme+1 ; j++){
		for (int k=0; k!=Nkterme+1; k++){
			tab[Niterme][j][k] = min(borne_fine, tab[Niterme][j][k]);
		}
	}
	
	for(int i = Niterme - 1; i != 0 ; i--){
		
		if(i%100==0){
			cout << "i= " << i  << endl;
		}

		for(int k =0; k!=Nkterme+1; k++){
			tab[i][Njterme][k]=tabOpt[Njterme][k];
		}

		for(int j =1;j!=Njterme+1; j++){
			tab[i][j][Nkterme]=tabOpt[j][Nkterme];
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

//tabtau4 est déjà optimisé
vector< vector< vector< vector<double> > > > tabTau4(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d ){

	cout << "debut"  << endl;
	vector < vector< vector<double> > > tabOpt(tabTau3Opt(NiTau2, NjTau2, NiTau3,  NjTau3, NkTau3, d));
	cout << "fin du calcul de tabTau3"  << endl;

	vector< vector < vector< vector<double> > > > tab(Niterme+1, vector< vector< vector<double> > >(Njterme+1, vector < vector<double> > (Nkterme+1, vector<double>(Nlterme+1))));

	double borne_fine;
	if(Niterme==1000){
		borne_fine=0.39;
	}else{
		borne_fine=0.20979;
	}

	//fin de l'initialisation du cube Nterme+1

	for(int j=Njterme; j!=0; j--){
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				tab[Niterme][j][k][l]=min(tabOpt[j][k][l], borne_fine);
			}
		}
	}

	tab[Niterme][0][0][0]=borne_fine;
	
	cout << "fin initialisation"  << endl;

	for(int i = Niterme - 1; i != 0 ; i--){
		//on print l'étape en cours
		if(i%10==0){
			cout << "i= " << i  << endl;
		}

		//On initialise chacune des faces 'loin' du cube i 
		for (int j=0; j!=Njterme+1; j++){
			for(int k=0; k!=Nkterme+1; k++){
				tab[i][j][k][Nlterme]=tabOpt[j][k][Nlterme];
			}
		}
		for (int k=0; k!=Nkterme+1; k++){
			for(int l=0; l!=Nlterme+1; l++){
				tab[i][Njterme][k][l]=tabOpt[Njterme][k][l];
			}
		}
		for (int j=0; j!=Njterme+1; j++){
			for(int l=0; l!=Nlterme+1; l++){
				tab[i][j][Nkterme][l]=tabOpt[j][Nkterme][l];
			}
		}
		//fin de l'initit=alisation

		for(int j = Njterme-1; j!=-1; j--){
			for(int k =Nkterme-1; k!=-1; k--){ 
				for (int l=Nlterme-1;l!=-1; l--){ 			
					double s_i = borneInfS(i, d);
					double s_j = borneInfS(j, d);
					double s_k = borneInfS(k, d);
					double s_l = borneInfS(l, d);
					double B_0_1 = max(i-j, 0);
					double B_1_2 = max(j-k, 0);
					double B_2_3 = max(k-l, 0);
					tab[i][j][k][l]=(1 + s_i * tab[i+1][j][k][l] + (s_j+B_0_1) * tab[i][j+1][k][l] + (s_k+B_1_2) * tab[i][j][k+1][l] + (s_l+B_2_3) * tab[i][j][k][l+1])/(s_i+s_j+s_k+s_l+B_0_1+B_1_2+B_2_3+l);
				}
			}
		}
    	}
	return tab;

}
double Tau4(int NiTau2, int NjTau2,int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d){
	vector< vector < vector< vector<double> > > > tab(tabTau4(NiTau2, NjTau2, NiTau3,NjTau3, NkTau3, Niterme, Njterme, Nkterme, Nlterme, d ));
	return tab[1][0][0][0];
}

vector< vector< vector< vector< vector<double> > > > > tabTau5(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d ){
	
	cout << "debut"  << endl;
	vector< vector < vector< vector<double> > > > tabOpt(tabTau4(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, d ));
	cout << "fin du calcul de tabTau4"  << endl;

	vector <vector< vector < vector< vector<double> > > > > tab(Niterme+1, vector< vector< vector< vector<double> > >>(Njterme+1, vector< vector < vector<double> > >(Nkterme+1, vector< vector<double> >(Nlterme+1, vector<double>(Nmterme +1)))));

	cout << "fin creation tabTau5"  << endl;
	double borne_fine;
	if(Niterme==100){
		borne_fine=1.17;
	}else if (Niterme ==1000){
		borne_fine=0.68;
	}else{
		borne_fine=0.04;
	}

	//fin de l'initialisation du cube Nterme+1
	for(int j=Njterme; j!=0; j--){
		cout << "j = " << j << endl;
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				for (int m=Nmterme; m!=-1; m--){
					tab[Niterme][j][k][l][m]=min(tabOpt[j][k][l][m], borne_fine);
				}
			}
		}
	}
	tab[Niterme][0][0][0][0]=borne_fine;
	cout << "fin initialisation tabTau5"  << endl;

	for(int i = Niterme - 1; i != 0 ; i--){
		
		if(i%10==0){
			cout << "i= " << i  << endl;//on print l'étape en cours
		}

		tab[i][Njterme]=tabOpt[Njterme];
		//cout << "fin initialisation cube i,j"  << endl;

		for(int j = Njterme-1; j!=-1; j--){
			//On initialise chacune des faces 'loin' du cube i, j
			//cout << j  << endl;

			for( int k=Nkterme; k!=-1; k--){
				for(int l=Nlterme; l!=-1; l--){
					//cout << "k = " << k << " l = " << l << endl;
					tab[i][j][k][l][Nmterme]=tab[Niterme][j][k][l][Nmterme];
				}
			}
			//cout << j  << endl;
			for(int k=Nkterme; k!=-1; k--){
				for(int m=Nmterme; m!=-1; m--){
					tab[i][j][k][Nlterme][m]=tab[Niterme][j][k][Nlterme][m];
				}
			}
			//cout << j  << endl;
			for(int l =Nlterme; l!=-1; l--){
				for(int m=Nmterme; m!=-1; m--){
					tab[i][j][Nkterme][l][m]=tab[Niterme][j][Nkterme][l][m];
				}
			}
			//cout << "fin initialisation cube k,i,j"  << endl;
			//fin de l'initialisation
			
			for(int k =Nkterme-1; k!=-1; k--){ 
				for (int l=Nlterme-1;l!=-1; l--){ 
					for(int m=Nmterme; m!=-1; m--){			
						double s_i = borneInfS(i, d);
						double s_j = borneInfS(j, d);
						double s_k = borneInfS(k, d);
						double s_l = borneInfS(l, d);
						double s_m = borneInfS(m, d);
						double B_0_1 = max(i-j, 0);
						double B_1_2 = max(j-k, 0);
						double B_2_3 = max(k-l, 0);
						double B_3_4 = max(l-m, 0);
						tab[i][j][k][l][m]=(1 + s_i * tab[i+1][j][k][l][m] + (s_j+B_0_1) * tab[i][j+1][k][l][m] + (s_k+B_1_2) * tab[i][j][k+1][l][m] + (s_l+B_2_3) * tab[i][j][k][l+1][m] + (s_m + B_3_4) * tab[i][j][k][l][m+1] )/(s_i + s_j + s_k + s_l + s_m + B_0_1 + B_1_2 + B_2_3 + B_3_4 + m);
					}
				}
			}
		}
    	}
	return tab;


}
double Tau5(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d ){

	vector <vector< vector < vector< vector<double> > > > > tab(tabTau5(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Niterme, Njterme, Nkterme, Nlterme, Nmterme, d));
	return tab[1][0][0][0][0];
}




