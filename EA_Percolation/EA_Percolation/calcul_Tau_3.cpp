#include <chrono>
#include "fonctions.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char * argv[]) {
	//int d = 30; //des choix atucieux de ni, nj, .. sont mis en commentaires
	int ni=100; //10000
	int nj=20; //1000
	int nk=30;//100 
	//int nl=10;
	int niTau2=1000;//10000
	int njTau2=100;//1000
	/*int niTau3 = 1000;
	int njTau3=1000;
	int nkTau3=100;*/
	int NiNeeded = 10;	
	int tabd[1] = {30};

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(int i=0; i<1; i++){
		int d=tabd[i];
	        
		//double res1 = Tau3Opt(niTau2, njTau2, ni, nj, nk, 1, d);
		//double res2 = Tau3Opt(niTau2, njTau2, ni, nj, nk, 5, d);
		//double res3 = Tau3Opt(niTau2, njTau2, ni, nj, nk, ni, d);
		//cout << "d = " << d  << " ni = " << ni << " nj = "<< nj <<" nk = "<< nk   << endl;
		//cout << "Tau_3 * d * 2 /log(2d) /3 = " << res * 2 * d /log(2*d) /3  << endl;
		cout << "Tau_3 / 3 = " << res1/3  << endl;
		cout << "borne diag = " << 0.3313/sqrt(d)  << endl;

		cout << "Tau3 = " << res1  << endl;
		//cout << "Tau3 = " << res2  << endl;
		cout << "Tau3 = " << res3  << endl;


	}
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
        cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
