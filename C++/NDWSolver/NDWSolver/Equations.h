#pragma once
#include<iostream>
#include<math.h>
#include<vector>
using namespace std;

class Equations {
private:
	double gamma;
	double a;
	double ndf1(double t, double state[]);
	double ndf2(double t, double state[]);
	bool abs;
public:
	int size();
	Equations(double g, double a0, bool absolute);
	double dx(double t, double state[]);
	double du(double t, double state[]);
	double dw(double t, double state[]);
};

Equations::Equations(double g=1., double a0=10., bool absolute = true)
{
	gamma = g;
	a = a0;
	abs = absolute;
};

double Equations::ndf1(double t, double state[])
{
	return 1. - exp(2. * state[1] - state[2]);
};

double Equations::ndf2(double t, double state[])
{
	return a * exp(-state[0] - state[2]) - 2.;
};

double Equations::dx(double t, double state[])
{
	if (abs) { return fabs(ndf1(t, state)); };
	return ndf1(t, state);
};

double Equations::du(double t, double state[])
{
	if (abs) { return fabs(ndf2(t, state)); };
	return ndf2(t, state);
};

double Equations::dw(double t, double state[])
{
	if (abs) {  
		return -(gamma - 1.) * (2. * fabs(ndf1(t, state)) + fabs(ndf2(t, state))); 
	};
	return -(gamma - 1.) * (2. * ndf1(t, state) + ndf2(t, state));
};

int Equations::size()
{
	return 3;
}