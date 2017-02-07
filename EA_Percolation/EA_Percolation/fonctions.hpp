//
//  fonctions.hpp
//  EA_Percolation
//
//  Created by Jules Pertinand on 07/02/2017.
//  Copyright © 2017 QBE & JPP. All rights reserved.
//

#ifndef fonctions_hpp
#define fonctions_hpp

#include <stdio.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

std::vector<double> tabTau1(int Nterme, int d);
double Tau1(int Nterme, int d);
double borneInfS(int i, int d);

std::vector< std::vector<double> > tabTau2(int Niterme, int Njterme, int d );
double Tau2(int Niterme, int Njterme, int d );

#endif /* fonctions_hpp */
