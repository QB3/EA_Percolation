//
//  fonctions.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include "fonctions.hpp"
using namespace std;
#include <math.h>  

//this function returns an array of size Nterme +1
//the i th element of the array is an bound on the quantity T_i
//T_0 is left empty (ie to 0)
//a rought bound on T_{Nterme} (Nterme \in \mathbb{N}) is computed
//then a tighter bound is obtained for T_{Nterme-1}, T_{Nterme-2}, ... and T_1 = \mathbb{E}(\tau_1)
vector<double> tabTau1(int Nterme, int dimension){
    vector<double> tab(Nterme+1); //the array of size Nterme +1 is created
    tab[Nterme] = 1.0/Nterme; // the array is initialized, ie T_{Nterme} is bound roughly by 1.0/Nterme
    for(int i = Nterme - 1; i != 0 ; i--){//the bound on T_i (stored in tab[i]) is computed backward with the formula presented in the article
        double s_i = borneInfS(i, dimension);//s_i is a bound on |S|
        tab[i] = (1.0 + tab[i+1]*s_i)/(s_i+i);//T_i is computed backward with the bound on T_{i+1}
    }
    return tab;//we return the whole tab
}

//this function encapsulates the previous one and returns the bound on T_1 (ie returns tab[1]) 
double Tau1(int Nterme, int d){
    vector<double> tab = tabTau1(Nterme, d);//call to the previous funciton
    return tab[1];//returns a bound on T_1
}

//this function computes the bound on |S| with the formula given by Dhar in his article (equation 7 in his paper)
double borneInfS(int i, int d){
	if(i==0){
		return 0;
	}
	if(d==2){
		return 2;
	}
    return ceil(2*(d-1)*pow(i, (d-2)/(double)(d-1)));
}

//returns an array of size Niterme+1 by Njterme +1
//the element [i, j] of the array is a bound on T_{i, j}
//T_{0,0} is left empty, ie left to 0
//The array T_{Niterme, j} is bound roughly by  T_{Niterme} 
//The array T_{i, Njterme} is bound roughly by  T_{Niterme}
//Then the Niterme+1 by Njterme +1 is filled with the recursiv inequality
vector< vector<double> > tabTau2(int Niterme, int Njterme, int d ){
    vector< vector<double> > tab(Niterme+1, vector<double>(Njterme+1));//creation of the Niterme+1 by Njterme +1 array
    tab[Niterme] = tabTau1(Njterme, d); //one border is bound
    double borne_fine=borne_Tau_i_0_0(Niterme);
    tab[Niterme][0]=borne_fine;

    //initialization more precise
	for(int j = 1; j != Njterme ; j++){
		tab[Niterme][j]=min(borne_fine, tab[Niterme][j]);
	}
    for(int i = Niterme - 1; i != 0 ; i--){
		tab[i][Njterme]=1.0/Njterme;//the second border is bound roughly
		for(int j = Njterme-1; j!=-1; j--){
			double s_i = borneInfS(i, d);//intermediate quantity
			double s_j = borneInfS(j, d);//intermediate quantity
			double B = max(i-j, 0);//intermediate quantity
			tab[i][j]=(1+s_i * tab[i+1][j] + (s_j + B)*tab[i][j+1])/(s_i+s_j+B+j);//we apply the backward inequality
		}
	}
        return tab;
}

//this function is the same as tabTau2 BUT:
//it computes the array starting from T_{Niterme} but does not store useless information
//example: while computing T_{1, 0, 0} at the next step you will have to bound the limits T_{i,j} in order to apply the backward inequatlity
//We compute T_{i, j} using an Niterme + 1* Njterme + 1  array BUT we only store an Njterme + 1 * NiNeeded +1 array
vector< vector<double> > tabTau2(int Niterme, int Njterme, int NiNeeded, int d ){

	//first part of the function, we compute backward the T_{i,j} without storing the results
    vector<double> tab_int(tabTau1(Njterme, d));
    double A=0.787064;
    double borne_fine=sqrt(3 *M_PI)/sqrt(Niterme)+3.0/(Niterme*A)*exp(-Niterme*A*A/3);
    tab_int[0]=borne_fine;
    //rafinement de l'initialisation
	for(int j = 1; j != Njterme ; j++){
		tab_int[j]=min(borne_fine, tab_int[j]);
	}
    for(int i = Niterme - 1; i != NiNeeded-1 ; i--){
		tab_int[Njterme]=1.0/Njterme;
		for(int j = Njterme-1; j!=-1; j--){
			double s_i = borneInfS(i, d);
			double s_j = borneInfS(j, d);
			double B = max(i-j, 0);
			double M_i=2*(d-2)*i+2;
			tab_int[j]=max((1+M_i * tab_int[j] + (s_j + B)*tab_int[j+1])/(M_i+s_j+B+j), (1+s_i * tab_int[j] + (s_j + B)*tab_int[j+1])/(s_i+s_j+B+j));
		}
	}

	//second part of the function, we compute backward the T_{i,j} AND we store the results
	vector< vector<double> > tabNeeded(NiNeeded+1, vector<double>(Njterme+1));
	tabNeeded[NiNeeded]=tab_int;
        for(int i = NiNeeded - 1; i != 0 ; i--){
		tabNeeded[i][Njterme]=1.0/Njterme;
		for(int j = Njterme-1; j!=-1; j--){
			double s_i = borneInfS(i, d);
			double s_j = borneInfS(j, d);
			double B = max(i-j, 0);
			double M_i=2*(d-2)*i+2;
			if(i==0){
				M_i=0;
			}
			double Q_s=(1+M_i * tabNeeded[i+1][j] + (s_j + B)*tabNeeded[i][j+1])/(M_i+s_j+B+j);
			double Q_M= (1+s_i * tabNeeded[i+1][j] + (s_j + B)*tabNeeded[i][j+1])/(s_i+s_j+B+j);
			tabNeeded[i][j]=max(Q_s, Q_M);
		}
	}
        return tabNeeded;
}

//this two functions encapsulates the previous one, without et with the optimization
double Tau2(int Niterme, int Njterme, int d ){
	vector< vector<double> > tab(tabTau2(Niterme, Njterme,2, d ));
	return tab[1][0];
}
double Tau2(int Niterme, int Njterme, int NiNeeded, int d ){
	vector< vector<double> > tab(tabTau2(Niterme, Njterme, NiNeeded,  d ));
	return tab[1][0];
}
^
//auxiliar function whixh cpy the array
vector< vector<double> > recopieTau2(int NiTau2, int NjTau2, int Njterme, int Nkterme, int d){
	vector< vector<double> > tabGrand(tabTau2(NiTau2, NjTau2, Njterme, d));

	vector< vector<double> >  tabPetit(Njterme+1, vector<double>(Nkterme+1));
	for( int j = 1; j!= Njterme+1; j++){
		for(int k = 0; k!=Nkterme+1; k++){
			tabPetit[j][k]=tabGrand[j][k];
		}
	}
	return tabPetit;
}

//returns a 3-dimensional array of size Niterme +1 * Njterme +1 * Nkterme + 1
//element [i,j,k] is a bound on T_{i, j, k}
//as the previous functions, the border are bounded roughly in order to be able to apply the backward inequality
vector< vector< vector<double> > > tabTau3(int Niterme, int Njterme, int Nkterme,
 int d ){
	vector < vector< vector<double> > > tab(Niterme+1, vector< vector<double> >(Njterme+1, vector<double>(Nkterme+1)));// creates a 2-DIMENSIONAL ARRAY
	tab[Niterme] = tabTau2(Njterme,Nkterme, d);
	//à modifier
	double A=0.39255;
	double borne_fine=2*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/8)/(Niterme* A*A/8);
	//double borne_fine=pow(12, 1.0/3)*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/12)/(Niterme *A*A/12);
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

//returns a 3-dimensional array of size Niterme +1 * Njterme +1 * Nkterme + 1
//element [i,j,k] is a bound on T_{i, j, k}
//as the previous functions, the border are bounded roughly in order to be able to apply the backward inequality
vector< vector< vector<double> > > tabTau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d ){
	vector< vector<double> > tabOpt(recopieTau2(NiTau2, NjTau2, Njterme, Nkterme, d));//cretas the array

	double A=0.39255;
	double borne_fine=2*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/8)/(Niterme* A*A/8);//computes a tight bound on T_{Niterme, 0, 0}
	tabOpt[0][0]=borne_fine;

	//first step of the function, computes the bounds without storing them all
	for(int i = Niterme - 1; i != NiNeeded-1 ; i--){	
		if(i%100==0){
			cout << "i= " << i  << endl;
		}
		double s_i = borneInfS(i, d);
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=-1; k--){ 
				double s_k = borneInfS(k, d);
				double B_1_2 = max(j-k, 0);
				tabOpt[j][k] =1 + s_i*tabOpt[j][k] + (s_j+B_0_1)*tabOpt[j+1][k] + (s_k+B_1_2)*tabOpt[j][k+1];//applies the recursion equality
				tabOpt[j][k] = tabOpt[j][k]/(s_i + s_j + B_0_1 + s_k + B_1_2);//applies the recursion equality
			}
		}
   	}
   	//second step of the function, computes the bounds and stores them, into a 3-dimensional array
	vector < vector< vector<double> > > tab(NiNeeded+1, vector< vector<double> >(Njterme+1, vector<double>(Nkterme+1)));// creates a 3-DIMENSIONAL ARRAY
	tab[NiNeeded]=tabOpt;
	for(int i = NiNeeded - 1; i !=0 ; i--){
		if(i%100==0){
			cout << "i= " << i  << endl;
		}
		double s_i = borneInfS(i, d);
		for(int k =0; k!=Nkterme+1; k++){
			//tab[i][Njterme][k]=tabOpt[Njterme][k];
			tab[i][Njterme][k]=tab[NiNeeded][Njterme][k];
		}
		for(int j =0;j!=Njterme+1; j++){
			tab[i][j][Nkterme]=tab[NiNeeded][j][Nkterme];
		}
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=-1; k--){ 
				double s_k = borneInfS(k, d);
				double B_1_2 = max(j-k, 0);
				tab[i][j][k] = 1 + s_i * tab[i+1][j][k] + (s_j+B_0_1) * tab[i][j+1][k] + (s_k + B_1_2) * tab[i][j][k+1];//applies the recursion equality
				tab[i][j][k] = tab[i][j][k]/(s_i + s_j +B_0_1 +s_k +B_1_2);//applies the recursion equality
			}
		}
    }
	return tab;
}

//encapsulates the  tabTau3Opt function
double Tau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d ){
	vector < vector< vector<double> > > tab(tabTau3Opt(NiTau2, NjTau2, Niterme, Njterme, Nkterme,1,  d ));
	return tab[1][0][0];
}
//encapsulates the  tabTau3Opt function
double Tau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d ){
	vector < vector< vector<double> > > tab(tabTau3Opt(NiTau2, NjTau2, Niterme, Njterme, Nkterme,NiNeeded,  d ));
	return tab[1][0][0];
}

//auxiliar function use to copy an array
vector< vector< vector<double> > > recopieTau3(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Njterme, int Nkterme, int Nlterme, int d ){
	vector < vector< vector<double> > > tabGrand(tabTau3Opt(NiTau2, NjTau2, NiTau3,  NjTau3, NkTau3, Njterme, d));
	vector < vector< vector<double> > > tabPetit(Njterme+1, vector < vector<double> > (Nkterme+1, vector<double>(Nlterme+1)));
	for(int j=Njterme; j!=0; j--){
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				tabPetit[j][k][l]=tabGrand[j][k][l];
			}
		}
	}
	return tabPetit;
}


//returns a 4-dimensional array of size Niterme +1 * Njterme +1 * Nkterme + 1 * Nlterme + 1
//element [i,j,k,l] is a bound on T_{i, j, k, l}
//as the previous functions, the border are bounded roughly in order to be able to apply the backward inequality
vector< vector< vector< vector<double> > > > tabTau4(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int NiNeeded, int d ){
	cout << "debut"  << endl;
	vector < vector< vector<double> > > tabOpt(recopieTau3(NiTau2, NjTau2, NiTau3,  NjTau3, NkTau3, Njterme, Nkterme, Nlterme, d));
	cout << "fin du calcul de tabTau3"  << endl;
	double borne_fine;
	if(Niterme==100){
		borne_fine=0.75;
	}else if(Niterme==1000){
		borne_fine=0.39;
	}else{
		borne_fine=0.20979;//borne pour 10 000
	}
	tabOpt[0][0][0]=borne_fine;	
	cout << "fin initialisation"  << endl;

	for(int i = Niterme - 1; i != NiNeeded-1 ; i--){
		//on print l'étape en cours
		if(i%100==0){
			cout << "i= " << i  << endl;
		}
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
					tabOpt[j][k][l]=(1 + s_i * tabOpt[j][k][l] + (s_j+B_0_1) * tabOpt[j+1][k][l] + (s_k+B_1_2) * tabOpt[j][k+1][l] + (s_l+B_2_3) * tabOpt[j][k][l+1])/(s_i+s_j+s_k+s_l+B_0_1+B_1_2+B_2_3+l);
				}
			}
		}
    	}

	vector< vector < vector< vector<double> > > > tab(NiNeeded+1, vector< vector< vector<double> > >(Njterme+1, vector < vector<double> > (Nkterme+1, vector<double>(Nlterme+1))));
	tab[NiNeeded]=tabOpt;
	//fin de l'initialisation du cube Nterme+1

	/*for(int j=Njterme; j!=0; j--){
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				tab[Niterme][j][k][l]=min(tabOpt[j][k][l], borne_fine);
			}
		}
	}*/
	for(int i = NiNeeded - 1; i != 0 ; i--){
		//on print l'étape en cours
		if(i%100==0){
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

//encapsulates the previous funciton
double Tau4(int NiTau2, int NjTau2,int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d){
	vector< vector < vector< vector<double> > > > tab(tabTau4(NiTau2, NjTau2, NiTau3,NjTau3, NkTau3, Niterme, Njterme, Nkterme, Nlterme, 1, d ));
	return tab[1][0][0][0];
}

//auxiliar function which copies an array
vector< vector< vector< vector<double> > > > recopieTau4(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d){
	vector< vector < vector< vector<double> > > > tabGrand(tabTau4(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Njterme, d ));
	vector< vector < vector< vector<double> > > >  tabPetit(Njterme+1, vector< vector < vector<double> > >(Nkterme+1, vector< vector<double> >(Nlterme+1, vector<double>(Nmterme +1))));
	for(int j=Njterme; j!=0; j--){
		cout << "j = " << j << endl;
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				for (int m=Nmterme; m!=-1; m--){
					tabPetit[j][k][l][m]=tabGrand[j][k][l][m];
				}
			}
		}
	}
	return tabPetit;
}

//returns a 5-dimensional array of size Niterme +1 * Njterme +1 * Nkterme + 1 * Nlterme + 1 * Nmterme + 1
//element [i,j,k,l,m] is a bound on T_{i, j, k, l,m}
//as the previous functions, the border are bounded roughly in order to be able to apply the backward inequality
vector< vector< vector< vector< vector<double> > > > > tabTau5(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int NiNeeded, int d ){
	
	cout << "debut"  << endl;
	vector< vector < vector< vector<double> > > > tabOpt(recopieTau4(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Njterme, Nkterme, Nlterme, Nmterme, d ));
	cout << "fin du calcul de tabTau4"  << endl;

	double borne_fine;
	if(Niterme==100){
		borne_fine=1.17;
	}else if (Niterme ==1000){
		borne_fine=0.68;
	}else{
		borne_fine=0.4;//borne pour 10 000 //tester 100 000
	}
	tabOpt[0][0][0][0]=borne_fine;

	for(int i = Niterme - 1; i != NiNeeded-1 ; i--){
		
		if(i%10==0){
			cout << "i= " << i  << endl;//on print l'étape en cours
		}

		for(int j = Njterme-1; j!=-1; j--){
			
			for(int k =Nkterme-1; k!=-1; k--){ 
				for (int l=Nlterme-1;l!=-1; l--){ 
					for(int m=Nmterme-1; m!=-1; m--){			
						double s_i = borneInfS(i, d);
						double s_j = borneInfS(j, d);
						double s_k = borneInfS(k, d);
						double s_l = borneInfS(l, d);
						double s_m = borneInfS(m, d);
						double B_0_1 = max(i-j, 0);
						double B_1_2 = max(j-k, 0);
						double B_2_3 = max(k-l, 0);
						double B_3_4 = max(l-m, 0);
						tabOpt[j][k][l][m]=(1 + s_i * tabOpt[j][k][l][m] + (s_j+B_0_1) * tabOpt[j+1][k][l][m] + (s_k+B_1_2) * tabOpt[j][k+1][l][m] + (s_l+B_2_3) * tabOpt[j][k][l+1][m] + (s_m + B_3_4) * tabOpt[j][k][l][m+1] )/(s_i + s_j + s_k + s_l + s_m + B_0_1 + B_1_2 + B_2_3 + B_3_4 + m);
					}
				}
			}
		}
    }

	vector <vector< vector < vector< vector<double> > > > > tab(NiNeeded+1, vector< vector< vector< vector<double> > >>(Njterme+1, vector< vector < vector<double> > >(Nkterme+1, vector< vector<double> >(Nlterme+1, vector<double>(Nmterme +1)))));
	cout << "fin creation tabTau5"  << endl;
	tab[NiNeeded]=tabOpt;

	for(int i = NiNeeded - 1; i != 0 ; i--){
		
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
					tab[i][j][k][l][Nmterme]=tab[NiNeeded][j][k][l][Nmterme];
				}
			}
			//cout << j  << endl;
			for(int k=Nkterme; k!=-1; k--){
				for(int m=Nmterme; m!=-1; m--){
					tab[i][j][k][Nlterme][m]=tab[NiNeeded][j][k][Nlterme][m];
				}
			}
			//cout << j  << endl;
			for(int l =Nlterme; l!=-1; l--){
				for(int m=Nmterme; m!=-1; m--){
					tab[i][j][Nkterme][l][m]=tab[NiNeeded][j][Nkterme][l][m];
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

//encapsulates the previous function
double Tau5(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d ){

	vector <vector< vector < vector< vector<double> > > > > tab(tabTau5(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Niterme, Njterme, Nkterme, Nlterme, Nmterme,1, d));
	return tab[1][0][0][0][0];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////// From there the code is no longer relevant for the main part of our article                  ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////// It is for a related work using similar ideas and briefy mentionnes at the end of our article////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double borne_Tau_i_0(int Niterme){
	double A=0.787064;
    double borne_fine=sqrt(3 *M_PI)/sqrt(Niterme)+3.0/(Niterme*A)*exp(-Niterme*A*A/3);
    return borne_fine;
}

//renvoie i=un tableau de taille Niterme+1 par Njterme +1
vector< vector<double> > tabDeltaTau12(int Niterme, int Njterme, int d ){
    vector< vector<double> > tab(Niterme+1, vector<double>(Njterme+1));
    tab[Niterme] = tabTau1(Njterme, d);
    double borne_fine=borne_Tau_i_0(Niterme);
    tab[Niterme][0]=borne_fine;
    //rafinement de l'initialisation
	for(int j = 1; j != Njterme+1 ; j++){
		tab[Niterme][j]=min(borne_fine, tab[Niterme][j]);
		tab[Niterme][j]= min(tab[Niterme][j], 1.0/j);
	}

    for(int i = Niterme - 1; i != 0 ; i--){
    	double borne_Tau_i= borne_Tau_i_0(i);
		tab[i][Njterme]=1.0/Njterme;
		tab[i][Njterme]=min(borne_Tau_i, tab[i][Njterme]);
		double s_i = borneInfS(i, d);
			for(int j = Njterme-1; j!=0; j--){
				double s_j = borneInfS(j, d);
				double B = max(i-j, 0);
				tab[i][j] = (1+s_i * tab[i+1][j] + (s_j + B)*tab[i][j+1])/(s_i+s_j+B+j);
				tab[i][j] = min(tab[i][j], 1.0/j);
				tab[i][j] = min(tab[i][j], borne_Tau_i);
			}
			tab[i][0]=(s_i * tab[i+1][0] + i*tab[i][1])/(s_i+i);
			tab[i][0]=min(tab[i][0],borne_Tau_i);
	}
        return tab;
}
double DeltaTau12(int Niterme, int Njterme, int NiNeeded,  int d ){
	vector< vector<double> > tab(tabDeltaTau12(Niterme, Njterme, NiNeeded, d));
	return tab[1][0];
}
double DeltaTau12(int Niterme, int Njterme, int d ){
	vector< vector<double> > tab(tabDeltaTau12(Niterme, Njterme,  d ));
	return tab[1][0];
}

//renvoie i=un tableau de taille Niterme+1 par Njterme +1
double maxCorner(double constanteNumerateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double tab_i_plus_1_j, double tab_i_j_plus_1, int j){
	double resultat = maxCorner(constanteNumerateur, 0, borneInfSi,  borneSupSi, borneInfSjPlusB, tab_i_plus_1_j, tab_i_j_plus_1, j);
	return resultat;
}
double maxCorner(double constanteNumerateur, double constanteDenominateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double tab_i_plus_1_j, double tab_i_j_plus_1, int j){
	double numerateurA =constanteNumerateur;
	numerateurA = numerateurA + borneInfSi * tab_i_plus_1_j;
	numerateurA = numerateurA + borneInfSjPlusB * tab_i_j_plus_1;
	double denominateurA = constanteDenominateur + borneInfSi + borneInfSjPlusB + j;

	double numerateurB = constanteNumerateur;
	numerateurB = numerateurB + borneSupSi * tab_i_plus_1_j;
	numerateurB = numerateurB + borneInfSjPlusB * tab_i_j_plus_1;
	double denominateurB = constanteDenominateur + borneSupSi + borneInfSjPlusB + j;

	return max(numerateurA/denominateurA, numerateurB/denominateurB);
}
vector< vector<double> > tabDeltaTau12(int Niterme, int Njterme, int NiNeeded, int d ){

    vector<double> tab_int(tabTau1(Njterme, d));
    double A=0.787064;
    double borne_fine=sqrt(3 *M_PI)/sqrt(Niterme)+3.0/(Niterme*A)*exp(-Niterme*A*A/3);
    tab_int[0]=borne_fine;

    //rafinement de l'initialisation
	for(int j = 1; j != Njterme ; j++){
		tab_int[j]=min(borne_fine, tab_int[j]);
		tab_int[j] = min(tab_int[j], 1.0/j);

	}
        for(int i = Niterme - 1; i != NiNeeded-1 ; i--){
		tab_int[Njterme]=1.0/Njterme;
		double s_i = borneInfS(i, d);
		double borne_Tau_i= borne_Tau_i_0(i);
		for(int j = Njterme-1; j!=0; j--){
			double s_j = borneInfS(j, d);
			double B = max(i-j, 0);
			tab_int[j]= (1+s_i * tab_int[j] + (s_j + B)*tab_int[j+1])/(s_i+s_j+B+j);
			tab_int[j]= min(tab_int[j], 1.0/j);
			tab_int[j]=min(tab_int[j], borne_Tau_i);
			//tab_int[j]=maxCorner(1, s_i, M_i, s_j + B, tab_int[j], tab_int[j+1], j);
		}
		//tab_int[0]=max((M_i * tab_int[0] + i*tab_int[1])/(M_i+i), (s_i * tab_int[0] + i*tab_int[1])/(s_i+i));
		tab_int[0]=(s_i * tab_int[0] + i*tab_int[1])/(s_i+i);
		tab_int[0]=min(tab_int[0], borne_Tau_i);
	}

	vector< vector<double> > tabNeeded(NiNeeded+1, vector<double>(Njterme+1));
	tabNeeded[NiNeeded]=tab_int;

        for(int i = NiNeeded - 1; i != 0 ; i--){
		tabNeeded[i][Njterme]=1.0/Njterme;
		double s_i = borneInfS(i, d);
		double borne_Tau_i= borne_Tau_i_0(i);
		for(int j = Njterme-1; j!=0; j--){
			double s_j = borneInfS(j, d);
			double B = max(i-j, 0);
			tabNeeded[i][j]= (1+s_i * tabNeeded[i+1][j] + (s_j + B)* tabNeeded[i][j+1])/(s_i+s_j+B+j);
			tabNeeded[i][j]=min(tabNeeded[i][j], 1.0/j);
			tabNeeded[i][j]=min(tabNeeded[i][j], borne_Tau_i);
			//tabNeeded[i][j]=maxCorner(1, s_i, M_i, s_j + B, tabNeeded[i+1][j], tabNeeded[i][j+1], j);		
		}
		tabNeeded[i][0]= (s_i * tabNeeded[i+1][0] + i * tabNeeded[i][1])/(s_i+i);
		tabNeeded[i][0]=min(tabNeeded[i][0], borne_Tau_i);
		//tabNeeded[i][0]=maxCorner(0, s_i, M_i, i, tabNeeded[i+1][0], tabNeeded[i][1], 0);
	}

    return tabNeeded;
}
vector< vector<double> > recopieDeltaTau12(int NiTau2, int NjTau2, int Njterme, int Nkterme, int d){
	vector< vector<double> > tabGrand(tabDeltaTau12(NiTau2, NjTau2, Njterme, d));

	vector< vector<double> >  tabPetit(Njterme+1, vector<double>(Nkterme+1));
	for( int j = 1; j!= Njterme+1; j++){
		for(int k = 0; k!=Nkterme+1; k++){
			tabPetit[j][k]=tabGrand[j][k];
		}
	}
	return tabPetit;
}
double maxCorner(double constanteNumerateur, double constanteDenominateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double borneSupSjPlusB, double borneInfSkPlusB, double tab_i_plus_1_j_k, double tab_i_j_plus_1_k, double tab_i_j_k_plus_1,int k){
	double NewConstNum = constanteNumerateur + borneInfSi * tab_i_plus_1_j_k;
	double NewConstDenom = constanteDenominateur + borneInfSi;
	double resultat1 = maxCorner(NewConstNum , NewConstDenom, borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,tab_i_j_plus_1_k, tab_i_j_k_plus_1, k);

	NewConstNum = constanteNumerateur + borneSupSi * tab_i_plus_1_j_k;
	NewConstDenom = constanteDenominateur + borneSupSi;
	double resultat2 = maxCorner(NewConstNum , NewConstDenom, borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,tab_i_j_plus_1_k, tab_i_j_k_plus_1, k);
	double resultat = max(resultat1, resultat2);
	return resultat;
}
double maxCorner(double constanteNumerateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double borneSupSjPlusB, double borneInfSkPlusB, double tab_i_plus_1_j_k, double tab_i_j_plus_1_k, double tab_i_j_k_plus_1,int k){
	double constanteDenominateur = 0;
	double NewConstNum = constanteNumerateur + borneInfSi * tab_i_plus_1_j_k;
	double NewConstDenom = constanteDenominateur + borneInfSi;
	double resultat1 = maxCorner(NewConstNum , NewConstDenom, borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,tab_i_j_plus_1_k, tab_i_j_k_plus_1, k);

	NewConstNum = constanteNumerateur + borneSupSi * tab_i_plus_1_j_k;
	NewConstDenom = constanteDenominateur + borneSupSi;
	double resultat2 = maxCorner(NewConstNum , NewConstDenom, borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,tab_i_j_plus_1_k, tab_i_j_k_plus_1, k);
	double resultat = max(resultat1, resultat2);
	return resultat;
}

double borne_Tau_i_0_0(int Niterme){
	double A=0.39255;
	double borne_fine=2*0.89298/(pow(Niterme, 1.0/3)) + exp(- Niterme * A*A*A/8)/(Niterme* A*A/8);
	return borne_fine;
}
//renvoie un tableau de taille Niterme + 1 * Njterme + 1 * Nkterme +1 * Nlterme + 1 
//à enventuellement rajouter, prendre le min avec la borne fine à chaque étape 
vector< vector< vector<double> > > tabDeltaTau23(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d ){
	//vector< vector<double> > tabOpt(tabTau2(NiTau2, NjTau2, d));
	vector< vector<double> > tabOpt(recopieDeltaTau12(NiTau2, NjTau2, Njterme, Nkterme, d));
	double borne_fine = borne_Tau_i_0_0(Niterme);
	tabOpt[0][0]=borne_fine;
	tabOpt[0][0]=min(tabOpt[1][0], tabOpt[0][0]);

	double numerateurA;
	double denominateurA;
	cout<<"borne fine : " << borne_fine << endl;
	for(int i = Niterme - 1; i != NiNeeded-1 ; i--){	
		if(i%100==0){
			cout << "i= " << i  << endl;
		}
		double s_i = borneInfS(i, d);
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=0; k--){ 
				double s_k = borneInfS(k, d);
				double B_1_2 = max(j-k, 0);
				numerateurA = 1;
				numerateurA = numerateurA + s_i * tabOpt[j][k];
				numerateurA = numerateurA + (s_j+B_0_1) * tabOpt[j+1][k];
				numerateurA = numerateurA +  (s_k+B_1_2) * tabOpt[j][k+1];
				denominateurA = s_i+s_j+s_k+B_0_1+B_1_2+k;
				tabOpt[j][k]=numerateurA/denominateurA;
				//tabOpt[j][k] = maxCorner(1, s_i, M_i, s_j+B_0_1 , M_j + i, s_k + B_1_2, tabOpt[j][k], tabOpt[j+1][k], tabOpt[j][k+1], k);
			}
			//tabOpt[j][0] = maxCorner(0, s_i, M_i, s_j+B_0_1 , M_j + i, j, tabOpt[j][0], tabOpt[j+1][0], tabOpt[j][1], 0);
			numerateurA = 0;
			numerateurA = numerateurA + s_i * tabOpt[j][0];
			numerateurA = numerateurA + (s_j+B_0_1) * tabOpt[j+1][0];
			numerateurA = numerateurA +  j * tabOpt[j][1];
			denominateurA = s_i+s_j+B_0_1+j;
			tabOpt[j][0]=numerateurA/denominateurA;
		}
    }

	vector < vector< vector<double> > > tab(NiNeeded+1, vector< vector<double> >(Njterme+1, vector<double>(Nkterme+1)));
	tab[NiNeeded]=tabOpt;

	for(int i = NiNeeded - 1; i !=0 ; i--){
		
		if(i%100==0){
			cout << "i= " << i  << endl;
		}

		for(int k =0; k!=Nkterme+1; k++){
			tab[i][Njterme][k]=tabOpt[Njterme][k];
		}

		for(int j =0;j!=Njterme+1; j++){
			tab[i][j][Nkterme]=tabOpt[j][Nkterme];
		}

		double s_i = borneInfS(i, d);
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=0; k--){ 
				double s_k = borneInfS(k, d);
				double B_1_2 = max(j-k, 0);
				//tab[i][j][k]=maxCorner(1, s_i, M_i, s_j+B_0_1 , M_j + i, s_k + B_1_2, tab[i+1][j][k], tab[i][j+1][k], tab[i][j][k+1], k);
				numerateurA = 1;
				numerateurA = numerateurA + s_i * tab[i+1][j][k];
				numerateurA = numerateurA + (s_j+B_0_1) * tab[i][j+1][k];
				numerateurA = numerateurA +  (s_k+B_1_2) * tab[i][j][k+1];
				denominateurA = s_i+s_j+s_k+B_0_1+B_1_2+k;
				tab[i][j][k]=numerateurA/denominateurA;
			}
			numerateurA = 0;
			numerateurA = numerateurA + s_i * tab[i+1][j][0];
			numerateurA = numerateurA + (s_j+B_0_1) * tab[i][j+1][0];
			numerateurA = numerateurA +  j * tab[i][j][1];
			denominateurA = s_i+s_j+B_0_1+j;
			tab[i][j][0]=numerateurA/denominateurA;
			//tab[i][j][0] = maxCorner(0, s_i, M_i, s_j+B_0_1 , M_j + i, j, tab[i+1][j][0], tab[i][j+1][0], tab[i][j][1], 0);
		}
    }
	return tab;
}
double DeltaTau23(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d ){
	vector < vector< vector<double> > > tab(tabDeltaTau23(NiTau2, NjTau2, Niterme, Njterme, Nkterme,1,  d ));
	return tab[1][0][0];
}
double DeltaTau23(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d ){
	vector < vector< vector<double> > > tab(tabDeltaTau23(NiTau2, NjTau2, Niterme, Njterme, Nkterme, NiNeeded,  d ));
	return tab[1][0][0];
}
vector< vector< vector<double> > > recopieDeltaTau23(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Njterme, int Nkterme, int Nlterme, int d ){
	vector < vector< vector<double> > > tabGrand(tabDeltaTau23(NiTau2, NjTau2, NiTau3,  NjTau3, NkTau3, Njterme, d));
	//cout << "delta tau 23= " << tabGrand[1][0][0]  << endl;
	vector < vector< vector<double> > > tabPetit(Njterme+1, vector < vector<double> > (Nkterme+1, vector<double>(Nlterme+1)));
	for(int j=Njterme; j!=0; j--){
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				tabPetit[j][k][l]=tabGrand[j][k][l];
			}
		}
	}
	return tabPetit;
}
double maxCorner(double constanteNumerateur, double constanteDenominateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double borneSupSjPlusB, double borneInfSkPlusB, double borneSupSkPlusB, double borneInfSlPlusB, double tab_i_plus_1_j_k_l, double tab_i_j_plus_1_k_l, double tab_i_j_k_plus_1_l, double tab_i_j_k_l_plus_1,int l){
	double NewConstNum = constanteNumerateur + borneInfSi * tab_i_plus_1_j_k_l;
	double NewConstDenom = constanteDenominateur + borneInfSi;
	double resultat1 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, tab_i_j_plus_1_k_l, tab_i_j_k_plus_1_l,  tab_i_j_k_l_plus_1,l);

	NewConstNum = constanteNumerateur + borneSupSi * tab_i_plus_1_j_k_l;
	NewConstDenom = constanteDenominateur + borneSupSi;
	double resultat2 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, tab_i_j_plus_1_k_l, tab_i_j_k_plus_1_l,  tab_i_j_k_l_plus_1,l);
	double resultat = max(resultat1, resultat2);
	return resultat;
}
double maxCorner(double constanteNumerateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double borneSupSjPlusB, double borneInfSkPlusB, double borneSupSkPlusB, double borneInfSlPlusB, double tab_i_plus_1_j_k_l, double tab_i_j_plus_1_k_l, double tab_i_j_k_plus_1_l, double tab_i_j_k_l_plus_1,int l){
	double constanteDenominateur=0;
	double NewConstNum = constanteNumerateur + borneInfSi * tab_i_plus_1_j_k_l;
	double NewConstDenom = constanteDenominateur + borneInfSi;
	double resultat1 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, tab_i_j_plus_1_k_l, tab_i_j_k_plus_1_l,  tab_i_j_k_l_plus_1,l);

	NewConstNum = constanteNumerateur + borneSupSi * tab_i_plus_1_j_k_l;
	NewConstDenom = constanteDenominateur + borneSupSi;
	double resultat2 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, tab_i_j_plus_1_k_l, tab_i_j_k_plus_1_l,  tab_i_j_k_l_plus_1,l);
	double resultat = max(resultat1, resultat2);
	return resultat;
}
vector< vector< vector< vector<double> > > > tabDeltaTau34(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int NiNeeded, int d ){

	cout << "debut"  << endl;
	vector < vector< vector<double> > > tabOpt(recopieDeltaTau23(NiTau2, NjTau2, NiTau3,  NjTau3, NkTau3, Njterme, Nkterme, Nlterme, d));
	//vector < vector< vector<double> > > tabOpt(recopieTau3(NiTau2, NjTau2, NiTau3,  NjTau3, NkTau3, Njterme, Nkterme, Nlterme, d));
	cout << "fin du calcul de tabTau3"  << endl;
	cout << "fin du calcul de tabTau3"  << endl;
	double borne_fine;
	if(Niterme==100){
		borne_fine=0.75;
	}else if(Niterme==1000){
		borne_fine=0.39;
	}else{
		borne_fine=0.20979;//borne pour 10 000
	}
	borne_fine=tabOpt[1][0][0];
	tabOpt[0][0][0]=borne_fine;	
	cout << "fin initialisation"  << endl;
	double numerateurA;
	double denominateurA;
	double numerateurB;
	double denominateurB;
	for(int i = Niterme - 1; i != NiNeeded-1 ; i--){
		//on print l'étape en cours
		if(i%100==0){
			cout << "i= " << i  << endl;
		}
		double s_i = borneInfS(i, d);
		double M_i = 2*(d-2)*i+2;
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double M_j = 2*(d-2)*j+2;
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=-1; k--){ 
				double s_k = borneInfS(k, d);
				double M_k = 2*(d-2)*k+2;
				double B_1_2 = max(j-k, 0);
				for (int l=Nlterme-1;l!=0; l--){ 			
					double s_l = borneInfS(l, d);
					double B_2_3 = max(k-l, 0);
					numerateurA = 1;
					numerateurA = numerateurA + s_i * tabOpt[j][k][l];
					numerateurA = numerateurA + (s_j+B_0_1) * tabOpt[j+1][k][l];
					numerateurA = numerateurA +  (s_k+B_1_2) * tabOpt[j][k+1][l];
					numerateurA = numerateurA +  (s_l+B_2_3) * tabOpt[j][k][l+1];
					denominateurA = s_i+s_j+s_k+s_l+B_0_1+B_1_2+B_2_3+l;
					tabOpt[j][k][l]=numerateurA/denominateurA;
				}
				numerateurA = 0;
				numerateurA = numerateurA + s_i * tabOpt[j][k][0];
				numerateurA = numerateurA + (s_j+B_0_1) * tabOpt[j+1][k][0];
				numerateurA = numerateurA +  (s_k+B_1_2) * tabOpt[j][k+1][0];
				numerateurA = numerateurA +  k * tabOpt[j][k][1];
				denominateurA = s_i+s_j+s_k+B_0_1+B_1_2+k;
				tabOpt[j][k][0]=numerateurA/denominateurA;
			}
		}
    }

	vector< vector < vector< vector<double> > > > tab(NiNeeded+1, vector< vector< vector<double> > >(Njterme+1, vector < vector<double> > (Nkterme+1, vector<double>(Nlterme+1))));
	tab[NiNeeded]=tabOpt;
	//fin de l'initialisation du cube Nterme+1
	for(int i = NiNeeded - 1; i != 0 ; i--){
		//on print l'étape en cours
		if(i%100==0){
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
		double s_i = borneInfS(i, d);
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double M_j = 2*(d-2)*j+2;
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=-1; k--){ 
				double s_k = borneInfS(k, d);
				double M_k = 2*(d-2)*k+2;
				double B_1_2 = max(j-k, 0);
				for (int l=Nlterme-1;l!=0; l--){ 			
					double s_l = borneInfS(l, d);
					double B_2_3 = max(k-l, 0);
					numerateurA = 1;
					numerateurA = numerateurA + s_i * tab[i+1][j][k][l];
					numerateurA = numerateurA + (s_j+B_0_1) * tab[i][j+1][k][l];
					numerateurA = numerateurA +  (s_k+B_1_2) * tab[i][j][k+1][l];
					numerateurA = numerateurA +  (s_l+B_2_3) * tab[i][j][k][l+1];
					denominateurA = s_i+s_j+s_k+s_l+B_0_1+B_1_2+B_2_3+l;
					tab[i][j][k][l]=numerateurA/denominateurA;
				}
				numerateurA = 0;
				numerateurA = numerateurA + s_i * tab[i+1][j][k][0];
				numerateurA = numerateurA + (s_j+B_0_1) * tab[i][j+1][k][0];
				numerateurA = numerateurA +  (s_k+B_1_2) * tab[i][j][k+1][0];
				numerateurA = numerateurA +  k * tab[i][j][k][1];
				denominateurA = s_i+s_j+s_k+B_0_1+B_1_2+k;
				tab[i][j][k][0] = numerateurA/denominateurA;
			}
		}
    }
	return tab;
}
double DeltaTau34(int NiTau2, int NjTau2,int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d){
	vector< vector < vector< vector<double> > > > tab(tabDeltaTau34(NiTau2, NjTau2, NiTau3,NjTau3, NkTau3, Niterme, Njterme, Nkterme, Nlterme, 1, d ));
	return tab[1][0][0][0];
}
double DeltaTau34(int NiTau2, int NjTau2,int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int NiNeeded, int d){
	vector< vector < vector< vector<double> > > > tab(tabDeltaTau34(NiTau2, NjTau2, NiTau3,NjTau3, NkTau3, Niterme, Njterme, Nkterme, Nlterme, 1, d ));
	return tab[1][0][0][0];
}
void print(vector< vector< double > > tab, int Njterme, int Nkterme){
	cout.precision(2);
	for(int j =0; j!= Njterme +1 ; j++){
		for(int k=0; k!=Nkterme+1; k++){
			cout << tab[j][k] << " ";	
		}
		cout<< endl;
	}
	cout<< endl;
}
vector< vector< vector< vector<double> > > > recopieDeltaTau34(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d){
	vector< vector < vector< vector<double> > > > tabGrand(tabDeltaTau34(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Njterme, d ));
	vector< vector < vector< vector<double> > > >  tabPetit(Njterme+1, vector< vector < vector<double> > >(Nkterme+1, vector< vector<double> >(Nlterme+1, vector<double>(Nmterme +1))));
	for(int j=Njterme; j!=0; j--){
		cout << "j = " << j << endl;
		for(int k = Nkterme; k!=-1; k--){
			for(int l =Nlterme; l!=-1; l--){
				for (int m=Nmterme; m!=-1; m--){
					tabPetit[j][k][l][m]=tabGrand[j][k][l][m];
				}
			}
		}
	}
	return tabPetit;
}
double maxCorner(double constanteNumerateur, double constanteDenominateur, double borneInfSi, double borneSupSi, double borneInfSjPlusB, double borneSupSjPlusB, double borneInfSkPlusB, double borneSupSkPlusB, double borneInfSlPlusB, double borneSupSlPlusB, double borneInfSmPlusB, double tab_i_plus_1_j_k_l_m, double tab_i_j_plus_1_k_l_m, double tab_i_j_k_plus_1_l_m, double tab_i_j_k_l_plus_1_m,double tab_i_j_k_l_m_plus_1 ,int m){
	double NewConstNum = constanteNumerateur + borneInfSi * tab_i_plus_1_j_k_l_m;
	double NewConstDenom = constanteDenominateur + borneInfSi;
	double resultat1 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, borneSupSlPlusB, borneInfSmPlusB, tab_i_j_plus_1_k_l_m, tab_i_j_k_plus_1_l_m, tab_i_j_k_l_plus_1_m, tab_i_j_k_l_m_plus_1 , m);

	NewConstNum = constanteNumerateur + borneSupSi * tab_i_plus_1_j_k_l_m;
	NewConstDenom = constanteDenominateur + borneSupSi;
	double resultat2 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, borneSupSlPlusB, borneInfSmPlusB, tab_i_j_plus_1_k_l_m, tab_i_j_k_plus_1_l_m, tab_i_j_k_l_plus_1_m, tab_i_j_k_l_m_plus_1 , m);
	double resultat = max(resultat1, resultat2);
	return resultat;
}
double maxCorner(double constanteNumerateur,  double borneInfSi, double borneSupSi, double borneInfSjPlusB, double borneSupSjPlusB, double borneInfSkPlusB, double borneSupSkPlusB, double borneInfSlPlusB, double borneSupSlPlusB, double borneInfSmPlusB, double tab_i_plus_1_j_k_l_m, double tab_i_j_plus_1_k_l_m, double tab_i_j_k_plus_1_l_m, double tab_i_j_k_l_plus_1_m,double tab_i_j_k_l_m_plus_1 ,int m){
	double constanteDenominateur = 0;
	double NewConstNum = constanteNumerateur + borneInfSi * tab_i_plus_1_j_k_l_m;
	double NewConstDenom = constanteDenominateur + borneInfSi;
	double resultat1 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, borneSupSlPlusB, borneInfSmPlusB, tab_i_j_plus_1_k_l_m, tab_i_j_k_plus_1_l_m, tab_i_j_k_l_plus_1_m, tab_i_j_k_l_m_plus_1 , m);

	NewConstNum = constanteNumerateur + borneSupSi * tab_i_plus_1_j_k_l_m;
	NewConstDenom = constanteDenominateur + borneSupSi;
	double resultat2 = maxCorner(NewConstNum , NewConstDenom,  borneInfSjPlusB, borneSupSjPlusB, borneInfSkPlusB,  borneSupSkPlusB, borneInfSlPlusB, borneSupSlPlusB, borneInfSmPlusB, tab_i_j_plus_1_k_l_m, tab_i_j_k_plus_1_l_m, tab_i_j_k_l_plus_1_m, tab_i_j_k_l_m_plus_1 , m);
	double resultat = max(resultat1, resultat2);
	return resultat;
}
vector< vector< vector< vector< vector<double> > > > > tabDeltaTau45(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int NiNeeded, int d ){
	
	cout << "debut"  << endl;
	vector< vector < vector< vector<double> > > > tabOpt(recopieDeltaTau34(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Njterme, Nkterme, Nlterme, Nmterme, d ));
	cout << "fin du calcul de tabTau4"  << endl;

	double borne_fine;
	if(Niterme==100){
		borne_fine=1.17;
	}else if (Niterme ==1000){
		borne_fine=0.68;
	}else{
		borne_fine=0.4;//borne pour 10 000 //tester 100 000
	}
	borne_fine=tabOpt[1][0][0][0];
	tabOpt[0][0][0][0]=borne_fine;


	double numerateurA;
	double denominateurA;
	cout  << "BORNE FINE = "<< borne_fine<<endl;
	for(int i = Niterme - 1; i != NiNeeded-1 ; i--){
		if(i%10==0){
			cout << "i= " << i  << endl;//on print l'étape en cours
		}
		double s_i = borneInfS(i, d);
		for(int j = Njterme-1; j!=-1; j--){
			double s_j = borneInfS(j, d);
			double B_0_1 = max(i-j, 0);
			for(int k =Nkterme-1; k!=-1; k--){
				double s_k = borneInfS(k, d);
				double B_1_2 = max(j-k, 0);
				for (int l=Nlterme-1;l!=-1; l--){
					double s_l = borneInfS(l, d);
					double B_2_3 = max(k-l, 0);
					for(int m=Nmterme-1; m!=0; m--){			
						double s_m = borneInfS(m, d);
						double B_3_4 = max(l-m, 0);
						numerateurA= 1 ;
						numerateurA=numerateurA + s_i * tabOpt[j][k][l][m];
						numerateurA=numerateurA + (s_j+B_0_1) * tabOpt[j+1][k][l][m] ;
						numerateurA=numerateurA + (s_k+B_1_2) * tabOpt[j][k+1][l][m] ;
						numerateurA=numerateurA + (s_l+B_2_3) * tabOpt[j][k][l+1][m] ;
						numerateurA=numerateurA + (s_m + B_3_4) * tabOpt[j][k][l][m+1];
						tabOpt[j][k][l][m]=numerateurA/denominateurA;
					}
					numerateurA= 0 ;
					numerateurA=numerateurA + s_i * tabOpt[j][k][l][0];
					numerateurA=numerateurA + (s_j+B_0_1) * tabOpt[j+1][k][l][0] ;
					numerateurA=numerateurA + (s_k+B_1_2) * tabOpt[j][k+1][l][0] ;
					numerateurA=numerateurA + (s_l+B_2_3) * tabOpt[j][k][l+1][0] ;
					numerateurA=numerateurA + l * tabOpt[j][k][l][1];
					denominateurA=s_i + s_j + s_k + s_l + B_0_1 + B_1_2 + B_2_3 + l;
					tabOpt[j][k][l][0]=numerateurA/denominateurA;
				}
			}
		}
   	}
	vector <vector< vector < vector< vector<double> > > > > tab(NiNeeded+1, vector< vector< vector< vector<double> > >>(Njterme+1, vector< vector < vector<double> > >(Nkterme+1, vector< vector<double> >(Nlterme+1, vector<double>(Nmterme +1)))));
	cout << "fin creation tabTau5"  << endl;
	tab[NiNeeded]=tabOpt;
	for(int i = NiNeeded - 1; i != 0 ; i--){
		double s_i = borneInfS(i, d);
		if(i%10==0){
			cout << "i= " << i  << endl;//on print l'étape en cours
		}
		tab[i][Njterme]=tabOpt[Njterme];
		for(int j = Njterme-1; j!=-1; j--){
			//On initialise chacune des faces 'loin' du cube i, j
			double s_j = borneInfS(j, d);
			double B_0_1 = max(i-j, 0);
			for( int k=Nkterme; k!=-1; k--){
				for(int l=Nlterme; l!=-1; l--){
					tab[i][j][k][l][Nmterme]=tab[NiNeeded][j][k][l][Nmterme];
				}
			}
			for(int k=Nkterme; k!=-1; k--){
				for(int m=Nmterme; m!=-1; m--){
					tab[i][j][k][Nlterme][m]=tab[NiNeeded][j][k][Nlterme][m];
				}
			}
			for(int l =Nlterme; l!=-1; l--){
				for(int m=Nmterme; m!=-1; m--){
					tab[i][j][Nkterme][l][m]=tab[NiNeeded][j][Nkterme][l][m];
				}
			}
			for(int k =Nkterme-1; k!=-1; k--){
				double s_k = borneInfS(k, d); 
				double M_k = 2*(d-2)*k+2 + j;
				double B_1_2 = max(j-k, 0);
				for (int l=Nlterme-1;l!=-1; l--){
					double s_l = borneInfS(l, d);
					double M_l = 2*(d-2)*l+2 + k;
					double B_2_3 = max(k-l, 0);
					for(int m=Nmterme-1; m!=0; m--){			
						double s_m = borneInfS(m, d);
						double B_3_4 = max(l-m, 0);
						numerateurA = 1 ;
						numerateurA = numerateurA + (s_i) * tab[i+1][j][k][l][m];
						numerateurA = numerateurA + (s_j+B_0_1) * tab[i][j+1][k][l][m] ;
						numerateurA = numerateurA + (s_k+B_1_2) * tab[i][j][k+1][l][m] ;
						numerateurA = numerateurA + (s_l+B_2_3) * tab[i][j][k][l+1][m] ;
						numerateurA = numerateurA + (s_m + B_3_4) * tab[i][j][k][l][m+1];
						denominateurA = s_i + s_j + s_k + s_l + s_m + B_0_1 + B_1_2 + B_2_3 + B_3_4 + m;
						tab[i][j][k][l][m]=numerateurA/denominateurA;
					}
					numerateurA = 0 ;
					numerateurA = numerateurA + s_i * tab[i+1][j][k][l][0];
					numerateurA = numerateurA + (s_j+B_0_1) * tab[i][j+1][k][l][0] ;
					numerateurA = numerateurA + (s_k+B_1_2) * tab[i][j][k+1][l][0] ;
					numerateurA = numerateurA + (s_l+B_2_3) * tab[i][j][k][l+1][0] ;
					numerateurA = numerateurA + l * tab[i][j][k][l][1];
					denominateurA = s_i + s_j + s_k + s_l + B_0_1 + B_1_2 + B_2_3 + l;
					tab[i][j][k][l][0]=numerateurA/denominateurA;
				}
			}
		}
    }
	return tab;
}
double DeltaTau45(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d ){
	vector <vector< vector < vector< vector<double> > > > > tab(tabDeltaTau45(NiTau2, NjTau2, NiTau3, NjTau3, NkTau3, NiTau4, NjTau4, NkTau4, NlTau4, Niterme, Njterme, Nkterme, Nlterme, Nmterme,1, d));
	return tab[1][0][0][0][0];
}