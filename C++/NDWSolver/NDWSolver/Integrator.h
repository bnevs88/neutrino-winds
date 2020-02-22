#pragma once
#include<iostream>
#include<math.h>
#include<vector>
#include "Equations.h"
using namespace std;

class Integrator {
private:
	vector<double> xsol; 
	vector<double> usol;
	vector<double> wsol;
	vector<double> tarray;
	Equations eqs;
	int size;
	
public:
	double* coupledRKStep(double t, double dt0, double state[]); //the output will be a pointer to an array
	Integrator(double initialState[], Equations equations);
	vector<double*> generateFunction(double initialState[], double xrange, double urange, int itermax); 
};

Integrator::Integrator(double initialState[], Equations equations)
{
	xsol = *new vector<double>(1, initialState[0]);
	usol = *new vector<double>(1, initialState[1]);
	wsol = *new vector<double>(1, initialState[2]);
	tarray = *new vector<double>(1, 0);
	eqs = equations;
}

double* Integrator::coupledRKStep(double t, double dt0, double state[])//adaptive RK4 step
{
	double karr[3][4] = { { 0,0,0,0 },{ 0,0,0,0 },{ 0,0,0,0 } };
	double* step=new double[3];
	double dt = dt0;
	double pc = 100;
	double statenorm = 0.;
	for (int i = 0; i < 3; i++) {
		statenorm += pow(state[i],2);
	}
	statenorm = sqrt(statenorm);
	while (pc > 1 or pc < .1) {
		for (int j = 0; j < 4; j++) { //j iterates through k1,k2,k3,k4 for the RK4 step
			if (j == 0) {
				karr[0][j] = dt * eqs.dx(t, state);
				karr[1][j] = dt * eqs.du(t, state);
				karr[2][j] = dt * eqs.dw(t, state);
			}
			else {
				//generate new state array to pass to functions
				double newstate[3] = { 0,0,0 };
				for (int k = 0; k < 3; k++) {
					newstate[k] = state[k] + karr[k][j - 1] / 2;
				}

				karr[0][j] = dt * eqs.dx(t + dt / 2, newstate);
				karr[1][j] = dt * eqs.du(t + dt / 2, newstate);
				karr[2][j] = dt * eqs.dw(t + dt / 2, newstate);
				for (int i = 0; i < 3; i++) {
					step[i] = state[i] + (karr[i][0] + 2 * karr[i][1] + 2 * karr[i][2] + karr[i][3]) / 6;
				}
			}
		}
		
		double stepnorm = 0.;
		for (int i = 0; i < 3; i++) {
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
	vector<double*> data;
	data.push_back(initialState);
	int i = 0;//iterator
	while (exp(data.back()[0]) > .0000001 && exp(data.back()[0]) < xrange && i < itermax) {
		//cout << "x: " << exp(data.back()[0]);
		data.push_back(coupledRKStep(0, .01, data.back()));
		i++;
	}
	return data;
}