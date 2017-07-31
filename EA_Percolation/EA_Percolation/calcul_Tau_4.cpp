//
//  main.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#include <chrono>
#include "fonctions.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char * argv[]) {
	//int d = 25;
	int ni=1000;//1000 à tester : 10000
	int nj=100;//100 à tester : 1000
	int nk=10;//100 à tester : 100
	int nl=10;//10
	int niTau2=10000;//10000
	int njTau2=1000;//1000
	int niTau3 = 1000;//1000
	int njTau3=100;//100
	int nkTau3=10;//10
	int tabd[8] = {25};
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
	for( int i =0; i< 1; i++){
		int d = tabd[i];
		double res = Tau4(niTau2, njTau2, niTau3, njTau3, nkTau3, ni, nj, nk, nl, d);
		cout << "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << " nl = " << nl << endl;
		cout << "Tau_4 * d * 2 /log(2d) /4 = " << res * d /log(2*d) /2  << endl;
		cout << "Tau_4 / 4 = " << res/4  << endl;
		cout << "borne diag = " << 0.3313/sqrt(d)  << endl;
	}
 
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
        cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
