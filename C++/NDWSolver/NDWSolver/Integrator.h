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
	double res;
	double vcrit;
	
public:
	double* coupledRKStep(double dt0, double state[], int itermax); //the output will be a pointer to an array
	Integrator(Equations equations, double resolution);
	vector<double*> generateFunction(double initialState[], double xrange, double urange, int itermax); 
	void writeToFile(string fileName);
	double findZeros(double initialState[], int itermax);
	double checkZeros(double initialState[], int itermax);
	double findV0(double tol, double maxprecision, int itermax);
	void scan(string fileName, double lower, double upper, int N, double tol);
	void printState(double state[]);
};

Integrator::Integrator(Equations equations, double resolution = 1.)
{
	eqs = equations;
	statelen = 4;
	res = resolution;
	vcrit = 0;
}

void Integrator::printState(double state[])
{
	for (int i = 0; i < 4; i++) {
		cout << state[i] << ", ";
	};
	cout << endl;
}

double* Integrator::coupledRKStep(double dt0, double state[], int itermax=100000)//adaptive RK4 step
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
	//cout << "Statenorm: " << statenorm << endl;
	int it = 0;
	while ((pc > res or pc < res/10.) and it<itermax) {
		for (int j = 0; j < 4; j++) { //j iterates through k1,k2,k3,k4 for the RK4 step
			if (j == 0) {
				karr[0][j] = dt * eqs.dx(state);
				karr[1][j] = dt * eqs.du(state);
				karr[2][j] = dt * eqs.dw(state);
			}
			else {
				//generate new state array to pass to functions
				double newstate[4] = { state[0]+dt/2.,0,0,0 };
				for (int k = 0; k < 3; k++) {//k iterates through x,u,w
					newstate[k+1] = state[k+1] + karr[k][j - 1] / 2;
				}

				karr[0][j] = dt * eqs.dx(newstate);
				karr[1][j] = dt * eqs.du(newstate);
				karr[2][j] = dt * eqs.dw(newstate);
				step[0] = state[0] + dt;
				//if (/*step[0] == .01875*/true) {
				//	cout << "karr" << endl;
				//	for (int i = 0; i < 3; i++) {
				//		for (int l = 0; l < 4; l++) {
				//			cout << karr[i][l] << endl;
				//		}
				//	}
				/*cout << "differentials: " << endl;
				cout << dt << endl;
				cout << eqs.dx(newstate) << endl;
				cout << eqs.du(newstate) << endl;
				cout << eqs.dw(newstate) << endl;
				cout << endl;*/
				//}
				for (int i = 0; i < 3; i++) {
					step[i+1] = state[i+1] + (karr[i][0] + 2. * karr[i][1] + 2. * karr[i][2] + karr[i][3]) / 6.;
				}
			}
		}
		
		double stepnorm = 0.;
		for (int i = 1; i < 4; i++) {
			stepnorm += pow(step[i] - state[i], 2);
		}
		stepnorm = sqrt(stepnorm);
		pc = 100 * stepnorm / statenorm;
		if (pc > res) { dt = dt / 2.; }
		else if (pc < res/10.) { dt = 1.9 * dt; };
		it++;
	}
	if (it >= itermax) { cout << "Max iteration count exceeded in coupledRKStep, pc=" << pc << endl; printState(step); };
	/*cout << "pc=" << pc << ", ";
	for (int i = 0; i < 4; i++) {
		cout << step[i] << ", ";
	};*/
	return step;
}

vector<double*> Integrator::generateFunction(double initialState[], double xrange=10, double urange=5, int itermax=100000)
{
	data.clear();
	data.push_back(initialState);
	int i = 0;//iterator
	//cout << "generating function" << endl;
	double* step = initialState;
	while (exp(step[1]) > .0000001 && exp(step[1]) < xrange && i < itermax) {
		//cout << "r=" << exp(data.back()[1]) << endl;
		step = coupledRKStep(.0001, data.back());
		data.push_back(step);
		i++;
	}
	if (i >= itermax) { cout << "Exceeded max iteration count in generateFunction" << endl; };
	//if (exp(data.back()[1]) <= .0000001) { cout << "x too small: " << exp(data.back()[1]) << endl; };
	//if (exp(data.back()[1]) >= xrange) { 
	//	cout << "x too big: "; 
	//	printState(data.end()[-1]); 
	//	cout << "dx: " << eqs.dx(data.end()[-1]) << endl;
	//	cout << "du: " << eqs.du(data.end()[-1]) << endl;
	//	cout << "dw: " << eqs.dw(data.end()[-1]) << endl;
	//};
	return data;
}

void Integrator::writeToFile(string fileName) 
{
	ofstream output(fileName);
	for (int i = 0; i < data.size(); i++) {
		for (int j = 0; j < 4; j++) {
			output << data[i][j]<<" ";
		}
		output << "\n";
	}
	output.close();
}

double Integrator::findZeros(double initialState[], int itermax=100000)
{
	data.clear();
	data.push_back(initialState);
	int i = 0;
	double tu = 0;
	double tx = 0;
	//cout << "finding zeros" << endl;
	while ((tx == 0 && tu == 0) && exp(data.back()[1])<500) {
		//cout << "r=" << exp(data.back()[1]) << endl;
		double* step = coupledRKStep(.01, data.back());
		/*cout << "state for f1: ";
		printState(data.back());
		cout << "step for f1: ";
		printState(step);
		cout << "current f1: " << eqs.ndf1(data.back()) << ", next f1: " << eqs.ndf1(step) << endl;*/
		if (copysign(1., eqs.ndf1(step)) != copysign(1., eqs.ndf1(data.back()))) {
			tu = step[0];
			//cout << "found tu=" << tu << endl;
		};
		if (copysign(1., eqs.ndf2(step)) != copysign(1., eqs.ndf2(data.back()))) {
			tx = step[0];
			//cout << "found tx=" << tx << endl;
		};
		data.push_back(step);
		i++;
	};
	if (tx==0 && tu==0) { 
		cout << "Exceeded max integration length in findZeros" << endl; 
		/*for (int j = 0; j < 4; j++) {
			cout << initialState[j] << endl;
		}*/
	}
	//cout << "tx: " << tx << ", tu: " << tu << ", tu-tx=" << tu - tx << endl;
	return tu - tx;
}

double Integrator::checkZeros(double initialState[], int itermax = 100000) 
{
	data.clear();
	data.push_back(initialState);
	int i = 0;
	double tu = 0;
	double tx = 0;
	//cout << "finding zeros" << endl;
	while ((tx == 0 || tu == 0) && exp(data.back()[1]) < 500) {
		//cout << "r=" << exp(data.back()[1]) << endl;
		double* step = coupledRKStep(.01, data.back());
		/*cout << "state for f1: ";
		printState(data.back());
		cout << "step for f1: ";
		printState(step);
		cout << "current f1: " << eqs.ndf1(data.back()) << ", next f1: " << eqs.ndf1(step) << endl;*/
		if (copysign(1., eqs.ndf1(step)) != copysign(1., eqs.ndf1(data.back()))) {
			tu = step[0];
			//cout << "found tu=" << tu << endl;
		};
		if (copysign(1., eqs.ndf2(step)) != copysign(1., eqs.ndf2(data.back()))) {
			tx = step[0];
			//cout << "found tx=" << tx << endl;
		};
		data.push_back(step);
		i++;
	};
	if (tx == 0 || tu == 0) {
		cout << "Exceeded max integration length in findZeros" << endl;
		/*for (int j = 0; j < 4; j++) {
			cout << initialState[j] << endl;
		}*/
	}
	//cout << "tx: " << tx << ", tu: " << tu << ", tu-tx=" << tu - tx << endl;
	return abs(tu - tx);
}

double Integrator::findV0(double tol = .9, double maxprecision = 1E-8, int itermax = 10000)
{
	if (vcrit == 0) {
		double v0 = .9;
		double dv = -.4;
		double testState[4] = { 0,0,0,0 };
		int i = 0;
		double err = 0;
		double prev = 0;
		bool exit = false;
		//cout << "Finding v0" << endl;
		do {
			while (abs(dv) > maxprecision&& i < itermax) {

				while ((copysign(1, findZeros(new double[4]{ 0,0,log(v0 + dv),0 })) != copysign(1, findZeros(new double[4]{ 0,0,log(v0),0 }))) or (v0 + dv > 1) or (v0 + dv < 0)) {
					//cout << "v0: " << v0 << ", dv: " << dv << endl;
					dv = dv / 2;
				};
				v0 = v0 + dv;
				//cout << "v0: " << v0 << ", dv: " << dv << endl;
				//cout << "tu-tx = " << findZeros(new double[4]{ 0, 0, log(v0), 0 }) << endl;
				i++;

			};
			err = checkZeros(new double[4]{ 0,0,log(v0 + dv),0 });
			if (err > tol) {
				if (100 * abs(prev - v0) / v0 > .5) {
					res = res / 2.;
					cout << "v0 is " << v0 + dv << ", error is " << err << ", increasing resolution to " << res << endl;
					prev = v0 + dv;
					v0 = .9;
					dv = -.4;
				}
				else {
					exit = true;
				};
			}
			else {
				exit = true;
			};
		} while (!exit);
		if (i >= itermax) { cout << "Max iteration count exceeded in findV0" << endl; };
		vcrit = v0 + dv;
	}
	return vcrit;
}

void Integrator::scan(string fileName, double lower = .001, double upper = .01, int N = 10, double tol = .9)
{
	ofstream output(fileName);
	for (double v = lower; v < upper; v = v + (upper-lower)/N) {
		cout << "Making data for v=" << v << endl;
		vector<double*> d = generateFunction(new double[4]{ 0,0,log(v),0 },50.);
		for (int i = 0; i < d.size(); i++) {
			for (int j = 0; j < 4; j++) {
				output << d[i][j] << " ";
			}
			output << "\n";
		}
		d.clear();
	}
	double v0 = findV0(tol);
	cout << "Making data for v=" << v0 << endl;
	vector<double*> d = generateFunction(new double[4]{ 0,0,log(v0),0 }, 50.);
	for (int i = 0; i < d.size(); i++) {
		for (int j = 0; j < 4; j++) {
			output << d[i][j] << " ";
		}
		output << "\n";
	}
	output.close();
}