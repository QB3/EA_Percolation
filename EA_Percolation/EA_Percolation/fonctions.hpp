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

using namespace std;

vector<double> tabTau1(int Nterme, int d);
double Tau1(int Nterme, int d);
double borneInfS(int i, int d);

vector< vector<double> > tabTau2(int Niterme, int Njterme, int d );
vector< vector<double> > tabTau2(int Niterme, int Njterme, int NiNeeded, int d );
double Tau2(int Niterme, int Njterme, int d );
double Tau2(int Niterme, int Njterme, int NiNeeded, int d );
vector< vector<double> > recopieTau2(int NiTau2, int NjTau2, int Njterme, int Nkterme, int d);

vector< vector< vector<double> > > tabTau3(int Niterme, int Njterme, int Nkterme, int d );
double Tau3(int Niterme, int Njterme, int Nkterme, int d );

vector< vector< vector<double> > > tabTau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d );
double Tau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d );
double Tau3Opt(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d );
vector< vector< vector<double> > > recopieTau3(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Njterme, int Nkterme, int Nlterme, int d );

vector< vector< vector< vector<double> > > > tabTau4(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d );
double Tau4(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d);
vector< vector< vector< vector<double> > > > recopieTau4(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Njterme, int Nkterme, int Nlterme, int Nmterme, int NiNeeded, int d);

vector< vector< vector< vector< vector<double> > > > > tabTau5(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d );
double Tau5(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d );

vector< vector<double> > tabDeltaTau12(int Niterme, int Njterme, int d );
vector< vector<double> > tabTau2(int Niterme, int Njterme, int NiNeeded, int d );
double DeltaTau12(int Niterme, int Njterme, int d );
double DeltaTau12(int Niterme, int Njterme, int NiNeeded,  int d );
vector< vector< vector<double> > > tabDeltaTau23(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d );
double DeltaTau23(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int d );
double DeltaTau23(int NiTau2, int NjTau2, int Niterme, int Njterme, int Nkterme, int NiNeeded, int d );
vector< vector<double> > recopieDeltaTau12(int NiTau2, int NjTau2, int Njterme, int Nkterme, int d);
vector< vector< vector< vector<double> > > > tabDeltaTau34(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int NiNeeded, int d );
double DeltaTau34(int NiTau2, int NjTau2,int NiTau3, int NjTau3, int NkTau3, int Niterme, int Njterme, int Nkterme, int Nlterme, int d);

vector< vector< vector< vector< vector<double> > > > > tabDeltaTau45(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int NiNeeded, int d );
vector< vector< vector< vector<double> > > > recopieDeltaTau34(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d);
double DeltaTau45(int NiTau2, int NjTau2, int NiTau3, int NjTau3, int NkTau3, int NiTau4, int NjTau4, int NkTau4, int NlTau4, int Niterme, int Njterme, int Nkterme, int Nlterme, int Nmterme, int d );

void print(vector< vector< double > > tab, int Njterme, int Nkterme);

#endif /* fonctions_hpp */
