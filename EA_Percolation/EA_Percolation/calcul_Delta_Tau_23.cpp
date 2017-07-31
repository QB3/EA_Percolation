#include <chrono>
#include "fonctions.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char * argv[]) {
	//int d = 30; //des choix atucieux de ni, nj, .. sont mis en commentaires
	int ni=1000; //10000
	int nj=100; //1000
	int nk=100;//100 
	//int nl=10;
	int niTau2=1000;//10000
	int njTau2=100;//1000
	/*int niTau3 = 1000;
	int njTau3=1000;
	int nkTau3=100;*/
	
	int tabd[1] = {18};

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(int i=0; i<1; i++){
		int d=tabd[i];
	    double res = DeltaTau23(niTau2, njTau2, ni, nj, nk, d);
		/*double res1 = DeltaTau23(niTau2, njTau2, ni, nj, nk,1, d);
		double res2 = DeltaTau23(niTau2, njTau2, ni, nj, nk,10, d);

		cout<< "DeltaTau23 = " << res1 << endl;
		cout<< "DeltaTau23 = " << res2 << endl;*/
		cout << "d = " << d  << " ni = " << ni << " nj = "<< nj <<" nk = "<< nk   << endl;
		cout << "DeltaTau23 = " << res << endl;
		cout << "borne diag = " << 0.3313/sqrt(d)  << endl;
	}
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
        cout << "Durée en milliseconde  = " << duration << endl;
}
