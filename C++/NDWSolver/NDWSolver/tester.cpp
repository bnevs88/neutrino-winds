#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<ctime>
#include "Equations.h"
#include "Integrator.h"
using namespace std;

int main() {
	cout << "\nHello world\n";
	Equations eqs(5./3.,.1);
	//Equations eqs;
	Integrator In(eqs);
	clock_t timereq = clock();
	double v0 = In.findV0();
	timereq = clock() - timereq;
	cout << "v0 = " << v0 << ", time elapsed = " << (float)timereq/CLOCKS_PER_SEC << endl;
	vector<double*> data = In.generateFunction(new double[4] { 0,0,log(v0),0 });
	cout << "Final state: " << endl;
	for (int i = 0; i < 4; i++) {
		cout << data.back()[i] << endl;
	}
	In.writeToFile("./output.txt");
	In.scan("./multiOutput.txt",.001,.005,10);
	
	return 0;
}

//trying to integrate without taking the absolute value results in an infinite loop that I haven't tracked down yet