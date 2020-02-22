#include<iostream>
#include<math.h>
#include<vector>
#include "Equations.h"
#include "Integrator.h"
using namespace std;

int main() {
	cout << "\nHello world\n";
	Equations eqs;
	double state[3] = { 0,log(.0050753194037),0 };
	Integrator In(state, eqs);
	cout << eqs.dx(0, state) <<", "<< eqs.du(0, state)<<", "<< eqs.dw(0, state);
	double *step = In.coupledRKStep(0, .1, state);
	cout << "\nStep:\n";
	for (int i = 0; i < 3; i++) {
		cout << step[i] << "\n";
	}

	double step2[3] = { 1.2,1.2,1.2 };
	double pc = 100;
	double statenorm = 0.;
	for (int i = 0; i < 3; i++) {
		statenorm += pow(state[i], 2);
	}
	statenorm = sqrt(statenorm);
	double stepnorm = 0.;
	
	for (int i = 0; i < 3; i++) {
		cout << "State: " << state[i] << endl;
		stepnorm += pow(step[i] - state[i], 2);
		cout << "added " << step[i] - state[i] << endl;
	}
	stepnorm = sqrt(stepnorm);
	pc = 100 * stepnorm / statenorm;
	cout << "statenorm: " << statenorm << endl;
	cout << "stepnorm: " << stepnorm << endl;
	cout << "Percent change: " << pc << endl;

	vector<double*> data = In.generateFunction(state);

	cout << "Final state: " << endl;
	for (int i = 0; i < 3; i++) {
		cout << data.back()[i] << endl;
	}
	return 0;
}