#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include "Equations.h"
#include "Integrator.h"
using namespace std;

int main() {
	cout << "\nHello world\n";
	Equations eqs;
	double state[4] = { 0,0,log(.01),0 };
	Integrator In(eqs);
	cout << eqs.dx(state) <<", "<< eqs.du(state)<<", "<< eqs.dw(state);
	double *step = In.coupledRKStep(.01, state);
	cout << "\nStep:\n";
	for (int i = 0; i < 4; i++) {
		cout << step[i] << "\n";
	}
	double v0 = In.findV0();
	cout << "v0 = " << v0 << endl;
	vector<double*> data = In.generateFunction(new double[4] { 0,0,log(v0),0 });

	cout << "Final state: " << endl;
	for (int i = 0; i < 4; i++) {
		cout << data.back()[i] << endl;
	}
	In.writeToFile("./output.txt");
	cout << "Zero at " << In.findZeros(state) << endl;
	
	return 0;
}