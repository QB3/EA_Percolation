//
//  main.cpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

//#include <chrono>
#include "fonctions.hpp"

using namespace std;
//using namespace std::chrono;

int main(int argc, const char * argv[]) {
    int d = 40;

    int ni=1000;
    int nj=1000;

   //return Tau1(n, d);
   // high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
	double res = Tau2(ni, nj, d);
     cout << "Tau_2 * d /log(2d) = " << res * d/log(2*d)  << endl;

    cout << "Tau_2 / 2 = " << res/2  << endl;

    cout << "borne diag = " << 0.3313/sqrt(d)  << endl;

   // high_resolution_clock::time_point t2 = high_resolution_clock::now();
    
    //auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    
    //cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
