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
    int d = 20;

    int ni=100;
    int nj=100;
   //return Tau1(n, d);
   // high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
    cout << "Tau 1 = " << Tau2(ni, nj, d) << endl;

   // high_resolution_clock::time_point t2 = high_resolution_clock::now();
    
    //auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    
    //cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
