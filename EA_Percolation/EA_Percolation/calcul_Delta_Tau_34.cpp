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
	int precision=4;
	int ni=1000;//1000 à tester : 10000
	int nj=100;//100 à tester : 1000
	int nk=10;//100 à tester : 100
	int nl=1;//10
	int niTau2=1000;//10000
	int njTau2=1000;//1000
	int niTau3 = 1000;//1000
	int njTau3=100;//100
	int nkTau3=100;//10
	//int tabd[10] = {2, 3, 4, 5, 10, 12, 15, 16, 17, 21};
	int tabd[1] = {16};


	int tailleTab=extent<decltype(tabd)>::value;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    fstream f("DeltaTau34.txt",ios_base::in | ios_base::out | ios_base::trunc);
	f.setf(std::ios::fixed,std::ios::floatfield);
	f.precision(precision);
	if (f.is_open()){
		for( int i =0; i< tailleTab; i++){
			int d = tabd[i];
			double res = DeltaTau34(niTau2, njTau2, niTau3, njTau3, nkTau3, ni, nj, nk, nl, d);
			cout << "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << " nl = " << nl << endl;
			cout.precision(4);
			cout << "DeltaTau34 = " << res  << endl;
			cout << "borne diag = " << 0.3313/sqrt(d)  << endl;
			f <<  d << "," <<  0.3313/sqrt(d)  << "," << res << endl;
	    }
	    f.close();
	}

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    cout << "Durée en milliseconde  = " << duration << endl;
}
