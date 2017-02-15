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
	int d = 22;
	int ni=1000;
	int nj=1000;//10-100
	int nk=10;//10-100
	int nl=10;//10-100
	int nm=10;
	int niTau2=10000;//10000
	int njTau2=1000;//1000
	int niTau3 = 1000;//1000-10000
	int njTau3=100;//1000
	int nkTau3=10;//100-1000
	int niTau4=1000;
	int njTau4=100;
	int nkTau4=10;
	int nlTau4=10;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
	double res = Tau5(niTau2, njTau2, niTau3, njTau3, nkTau3, niTau4, njTau4, nkTau4, nlTau4, ni, nj, nk, nl,nm, d);
	cout << "d = " << d  << " ni = " << ni << " nj = " << nj << " nk = "<< nk << " nl = " << nl << " nm = "<< nm << endl;
	cout << "Tau_5 * d * 2 /log(2d) /5 = " << res * 2 * d /log(2*d) /5  << endl;
	cout << "Tau_5 / 5 = " << res/5  << endl;
	cout << "borne diag = " << 0.3313/sqrt(d)  << endl;

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    
    cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
