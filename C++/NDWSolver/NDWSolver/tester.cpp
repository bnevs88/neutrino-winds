#include<iostream>
#include<math.h>
#include<vector>
#include "Equations.h"
#include "Integrator.h"
using namespace std;

int main() {
	cout << "\nHello world\n";
	Equations eqs;
	double state[3] = { 1,-2,1 };
	Integrator In(state, eqs);
	cout << eqs.dx(0, state) <<", "<< eqs.du(0, state)<<", "<< eqs.dw(0, state);
	double *step = In.coupledRKStep(0, .1, state);
	cout << "\nStep:\n";
	for (int i = 0; i < 3; i++) {
		cout << step[i] << "\n";
	}
	return 0;
}