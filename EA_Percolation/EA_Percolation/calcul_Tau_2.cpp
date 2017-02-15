#include <chrono>
#include "fonctions.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char * argv[]) {
	int d = 20;
	int ni=10000;
	int nj=1000;
	/*int nk=100;
	int nl=10;
	int niTau2=1000;
	int njTau2=1000;
	int niTau3 = 1000;
	int njTau3=1000;
	int nkTau3=100;*/

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    
	double res = Tau2(ni, nj, d);
	cout << "d = " << d  << " ni = " << ni << " nj = " << nj << endl;
	cout << "Tau_2 * d * 2 /log(2d) = " << res * 2 * d /log(2*d)  << endl;
	cout << "Tau_2 / 2 = " << res/2  << endl;

	cout << "borne diag = " << 0.3313/sqrt(d)  << endl;

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    
    auto duration = duration_cast<milliseconds>( t2 - t1 ).count();
    
    cout << "Durée en milliseconde  = " << duration << endl;
    
    //tabTau2(10, 10, 10);
    
    //return 0;
}
