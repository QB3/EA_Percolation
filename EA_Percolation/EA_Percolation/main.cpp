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
	int d = 40;
	int ni=1000;
	int nj=100;
	int nk=100;
	int niTau2=10000;
	int njTau2=1000;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
	double res = Tau3Opt(niTau2, njTau2, ni, nj, nk, d);
	cout << "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << endl;
	cout << "Tau_3 * d * 2 /log(2d) /3 = " << res * d * 2 /log(2*d) /3  << endl;
	cout << "Tau_3 / 3 = " << res/3  << endl;

	cout << "borne diag = " << 0.3313/sqrt(d)  << endl;

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    
    cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
