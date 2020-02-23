#pragma once
#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include "Equations.h"
using namespace std;

class Integrator {
private:
	vector<double*> data;
	Equations eqs;
	int statelen;
	
public:
	double* coupledRKStep(double dt0, double state[]); //the output will be a pointer to an array
	Integrator(Equations equations);
	vector<double*> generateFunction(double initialState[], double xrange, double urange, int itermax); 
	void writeToFile(string fileName);
	double findZeros(double initialState[], int itermax);
	double findV0(double maxprecision, int itermax);
};

Integrator::Integrator(Equations equations)
{
	eqs = equations;
	statelen = 4;
}

double* Integrator::coupledRKStep(double dt0, double state[])//adaptive RK4 step
{
	double karr[3][4] = { { 0,0,0,0 },{ 0,0,0,0 },{ 0,0,0,0 } };
	double* step=new double[4];
	double dt = dt0;
	double pc = 100;
	double statenorm = 0.;
	for (int i = 1; i < 4; i++) {
		statenorm += pow(state[i],2);
	}
	statenorm = sqrt(statenorm);
	while (pc > 1 or pc < .1) {
		for (int j = 0; j < 4; j++) { //j iterates through k1,k2,k3,k4 for the RK4 step
			if (j == 0) {
				karr[0][j] = dt * eqs.dx(state);
				karr[1][j] = dt * eqs.du(state);
				karr[2][j] = dt * eqs.dw(state);
			}
			else {
				//generate new state array to pass to functions
				double newstate[4] = { state[0]+dt/2.,0,0,0 };
				for (int k = 0; k < 3; k++) {
					newstate[k+1] = state[k+1] + karr[k][j - 1] / 2;
				}

				karr[0][j] = dt * eqs.dx(newstate);
				karr[1][j] = dt * eqs.du(newstate);
				karr[2][j] = dt * eqs.dw(newstate);
				step[0] = state[0] + dt;
				for (int i = 0; i < 3; i++) {
					step[i+1] = state[i+1] + (karr[i][0] + 2 * karr[i][1] + 2 * karr[i][2] + karr[i][3]) / 6;
				}
			}
		}
		
		double stepnorm = 0.;
		for (int i = 1; i < 4; i++) {
			stepnorm += pow(step[i] - state[i], 2);
		}
		stepnorm = sqrt(stepnorm);
		pc = 100 * stepnorm / statenorm;
		if (pc > 1) { dt = dt / 2; }
		else if (pc < .1) { dt = 2 * dt; };
	}
	/*cout << "\nStep:\n";
	for (int i = 0; i < 3; i++) {
		cout << step[i]<<"\n";
	}*/
	return step;
}

vector<double*> Integrator::generateFunction(double initialState[], double xrange=10, double urange=5, int itermax=10000)
{
	data.clear();
	data.push_back(initialState);
	int i = 0;//iterator
	while (exp(data.back()[1]) > .0000001 && exp(data.back()[1]) < xrange && i < itermax) {
		//cout << "x: " << exp(data.back()[0]);
		data.push_back(coupledRKStep(.01, data.back()));
		i++;
	}
	return data;
}

void Integrator::writeToFile(string filename) 
{
	ofstream output(filename);
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < 4; j++) {
			output << data[i][j]<<" ";
		}
		output << "\n";
	}
	output.close();
}

double Integrator::findZeros(double initialState[], int itermax=10000)
{
	data.clear();
	data.push_back(initialState);
	int i = 0;
	double tu = 0;
	double tx = 0;
	while (tx == 0 && tu == 0 && i < itermax) {
		double* step = coupledRKStep(.01, data.back());
		if (copysign(1., eqs.ndf1(step)) != copysign(1., eqs.ndf1(data.back()))) {
			tu = step[0];
		};
		if (copysign(1., eqs.ndf2(step)) != copysign(1., eqs.ndf2(data.back()))) {
			tx = step[0];
		};
		data.push_back(step);
	};
	return tu - tx;
}

double Integrator::findV0(double maxprecision = 1E-10, int itermax = 10000) {
	double v0 = 0;
	double dv = .5;
	double testState[4] = { 0,0,0,0 };
	int i = 0;
	while (dv > maxprecision&& i < itermax) {
		while (copysign(1, findZeros(new double[4] { 0,0,log(v0 + dv),0 }))!= copysign(1, findZeros(new double[4]{ 0,0,log(v0),0 }))) {
			dv = dv / 2;
		};
		v0 = v0 + dv;
		i++;
	};
	if (i >= itermax) { cout << "Max iteration count exceeded" << endl; };
	return v0 + dv;
}