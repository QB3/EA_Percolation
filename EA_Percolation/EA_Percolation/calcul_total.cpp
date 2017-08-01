//
//  main.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <chrono>
#include "fonctions.hpp"


using namespace std;
using namespace std::chrono;

int main(int argc, const char * argv[]) {
	//int d = 25;
	int precision=5;
	int ni=1000;
	int nj=100;//10-100
	int nk=10;//10-100
	int nl=10;//10-100
	int nm=1;
	int niTau2=1000;//10000
	int njTau2=1000;//1000
	int niTau3 = 1000;//1000-10000
	int njTau3=100;//1000
	int nkTau3=10;//100-1000
	int niTau4=1000;//10000
	int njTau4=100;//1000
	int nkTau4=10;//100
	int nlTau4=1;
	int tabd[11] = {2, 3, 4, 5, 10, 15, 20, 25, 30};
	//int tabd[4] = {2, 16, 18, 21};
	//int tabd[7] = {2, 16, 18, 21, 25, 30, 35};



	int tailleTab=extent<decltype(tabd)>::value;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    fstream f("calcul_total.txt",ios_base::in | ios_base::out | ios_base::trunc);
	f.setf(std::ios::fixed,std::ios::floatfield);
	f.precision(precision);
	if (f.is_open()){
		for( int i =0; i< tailleTab; i++){
			int d = tabd[i];
			double res2 = Tau2(niTau2, njTau2, d);
			double res3 = Tau3Opt(niTau2, njTau2, niTau3, njTau3, nkTau3,d);
			double res4 = Tau4(niTau2, njTau2, niTau3, njTau3, nkTau3, niTau4, njTau4, nkTau4, nlTau4, d);
			double res5 = Tau5(niTau2, njTau2, niTau3, njTau3, nkTau3, niTau4, njTau4, nkTau4, nlTau4, ni, nj, nk, nl,nm, d);
			cout << "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << " nl = " << nl << endl;
			cout.precision(4);
			//cout << "DeltaTau34 = " << res1  << endl;
			//cout << "borne diag = " << 0.3313/sqrt(d)  << endl;
			f <<  d << " & " <<  0.3313/sqrt(d)  << " & " << res2/2 << " & " << res3/3 << " & " << res4/4 << " & " << res5/5 <<  "\\" << "\\" << endl;
	    }
	    f.close();
	}

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    cout << "Durée en milliseconde  = " << duration << endl;
    //cout << "TAILLE TAB" << tailleTab;
}
