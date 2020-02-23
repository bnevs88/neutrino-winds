#pragma once
#include<iostream>
#include<math.h>
#include<vector>
using namespace std;

class Equations {
private:
	double gamma;
	double a;
	bool abs;
public:
	int size();
	Equations(double g, double a0, bool absolute);
	double ndf1(double state[]);
	double ndf2(double state[]);
	double dx(double state[]);
	double du(double state[]);
	double dw(double state[]);
};

Equations::Equations(double g=1., double a0=10., bool absolute = true)
{
	gamma = g;
	a = a0;
	abs = absolute;
};

double Equations::ndf1(double state[])
{
	return 1. - exp(2. * state[2] - state[3]);
};

double Equations::ndf2(double state[])
{
	return a * exp(-state[1] - state[3]) - 2.;
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
	if (abs) {  
		return -(gamma - 1.) * (2. * fabs(ndf1(state)) + fabs(ndf2(state))); 
	};
	return -(gamma - 1.) * (2. * ndf1(state) + ndf2(state));
};

int Equations::size()
{
	return 3;
}