#pragma once
#include<iostream>
#include<math.h>
#include<vector>
using namespace std;

class Equations {
private:
	double gamma;
	double a;
	double b;
	double betap;
	double lambda;
	bool abs;
public:
	int size();
	Equations(double g, double a0, double b0, double betap0, double lambda0, bool absolute);
	double ndf1(double state[]);
	double ndf2(double state[]);
	double dx(double state[]);
	double du(double state[]);
	double dw(double state[]);
};

Equations::Equations(double g = 1., double a0 = 10., double b0 = 0., double betap0 = 0., double lambda0=6., bool absolute = true)
{
	gamma = g;
	a = a0;
	b = b0;
	betap = betap0;
	lambda = lambda0;
	abs = absolute;
};

double Equations::ndf1(double state[])
{
	double f1 = 1. - exp(2. * state[2] - state[3]);
	if (isinf(f1)) { f1 = copysign(numeric_limits<double>::max(), f1); };
	return f1;
};

double Equations::ndf2(double state[])
{
	double f2 = a * exp(-state[1] - state[3]) - 2.;
	if (isinf(f2)) { f2 = copysign(numeric_limits<double>::max(), f2); };
	return f2;
};

double Equations::dx(double state[])
{
	if (abs) { return fabs(ndf1(state)); };
	return ndf1(state);
};

double Equations::du(double state[])
{
	if (abs) { return fabs(ndf2(state)); };
	return ndf2(state);
};

double Equations::dw(double state[])
{
	double f1 = ndf1(state);
	double f2 = ndf2(state);
	if (isinf(f1)) { f1 = copysign(numeric_limits<double>::max(), f1); };
	if (isinf(f2)) { f2 = copysign(numeric_limits<double>::max(), f2); };
	double dw1 = 0;
	if (abs) {
		dw1 = double(-(gamma - 1.) * (2. * fabs(f1) + fabs(f2) - b * fabs(f1) * exp(-state[1] - state[2] - state[3]) * (1 - betap * exp(2 * state[1] + lambda * state[3]))));
	}
	else {
		dw1 = double(-(gamma - 1.) * (2. * f1 + f2 - b * f1 * exp(-state[1] - state[2] - state[3]) * (1 - betap * exp(2 * state[1] + lambda * state[3]))));
	};
	if (!isnan(dw1)) { return dw1; }
	else {
		return copysign(numeric_limits<double>::max(), dw1);
	}; 
};

int Equations::size()
{
	return 3;
}