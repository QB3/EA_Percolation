//
//  main.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "fonctions.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char * argv[]) {
	//int d = 23;
	int ni=1000;
	int nj=100;//10-100
	int nk=10;//10-100
	int nl=10;//10-100
	int nm=10;
	int niTau2=1000;//10000
	int njTau2=1000;//1000
	int niTau3 = 1000;//1000-10000
	int njTau3=100;//1000
	int nkTau3=100;//100-1000
	int niTau4=1000;//10000
	int njTau4=100;//1000
	int nkTau4=100;//100
	int nlTau4=10;
	int tabd[8] = {2};
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
	for(int i =0; i<1; i++){
		int d=tabd[i];
		double res = DeltaTau45(niTau2, njTau2, niTau3, njTau3, nkTau3, niTau4, njTau4, nkTau4, nlTau4, ni, nj, nk, nl,nm, d);
		cout << "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << " nl = " << nl << " nm = "<< nm << endl;
		//cout << "Tau_5 * d * 2 /log(2d) /5 = " << res * 2 * d /log(2*d) /5  << endl;
		cout << "Delta_Tau_45 = " << res  << endl;
		cout << "borne diag = " << 0.3313/sqrt(d)  << endl;

		/*fstream f("tau5.txt",ios_base::in | ios_base::out | ios_base::trunc);
		if (f.is_open()){
			// Écrit les données :
			f<< "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << " nl = " << nl << " nm = "<< nm << endl << "Tau_5 / 5 = " << res/5  << endl<< "borne diag = " << 0.33133313/sqrt(d)  << endl;
			// Replace le pointeur de fichier au début :
			f.seekg(0);
			// Ferme le fichier :
			f.close();
	    	}*/

	}
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
        cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
