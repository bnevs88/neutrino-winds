#include <iostream> 
#include <fstream> 
#include <stdexcept>
#include <vector> 
#include <math.h>
#include <algorithm>
#include <ctime>
#include <cmath>

#include "integrator.h"
#include "rootfinder.h"
#include "EquationOfState.h"

const double PI = 2 * acos(0.0);

class GammaHeatingEquations : public EoS {
public:
	GammaHeatingEquations(double a, double b, double beta, double gamma, double v0, bool critStop = false) : a(a), b(b), beta(beta), gamma(gamma), v0(v0), rmax(100.), critStop(critStop) {}
	int size() { return 4; }
	std::vector<double> getInitialState() { return { 0.0,log(v0),0.0,0.0 }; }

	double f1(std::vector<double> state) { 
		double val = 1. - exp(2. * state[1])*pow(cs0,2)/pow(cs(state),2); 
		if (isinf(val)) { return copysign(LARGE_VALUE, val); };
		if (isnan(val)) { return copysign(LARGE_VALUE, val); };
		return val;
	};
	double f2(std::vector<double> state) { 
		//double val = a * exp(-state[0] - state[2]) - 2. + heat(state)*(gamma-1)/gamma;
		double val = a * exp(-state[0]) * pow(cs0, 2) / pow(cs(state), 2) - 2 + dpds(state) * heat(state) / (rho(state) * pow(cs(state), 2));
		if (isinf(val)) { return copysign(LARGE_VALUE, val); };
		if (isnan(val)) { return copysign(LARGE_VALUE, val); };
		return val;
	};
	double heat(std::vector<double> state) { 
		double val = b * exp(-state[0] - state[1] - state[2]) * (1.0 - beta * exp(2 * state[0] + 6.0 * state[2])); 
		if (isinf(val)) { return copysign(LARGE_VALUE, val); };
		if (isnan(val)) { return copysign(LARGE_VALUE, val); };
		return val;
	};

	std::vector<double> getDiffs(double t, std::vector<double> state) {
		double dw = -(gamma - 1) * (2. * f1(state) + f2(state) - f1(state) * heat(state));
		if (isinf(dw)) { dw = copysign(LARGE_VALUE, dw); };
		if (isnan(dw)) { dw = copysign(LARGE_VALUE, dw); };
		double dYe = 0.;
		//std::cout << "differentials: dx=" << f1(state) << ", du=" << f2(state) << ", dw=" << dw << std::endl;
		if (critStop) {
			return { f1(state), f2(state), dw, dYe };
		};
		return { fabs(f1(state)), fabs(f2(state)), dw, dYe };
	}

	double stopCondition(std::vector<double> state) {
		if (critStop) { return std::min(f1(state), f2(state)); };
		return rmax - exp(state[0]);
	}

	double ZeroFunction(double t, std::vector<double> state) {//return f1-f2, used to search for critical point
		return f1(state) - f2(state);
	}

protected:
	double a;
	double b;
	double beta;
	double gamma;
	double v0;
	double rmax;
	
	const double Mdot = 1.;
	const double r0 = 1.;
	const double cs0 = 1.;
	const double T0 = 1.;
	const double mb = 1.;

	bool critStop;
	const double LARGE_VALUE = 1.e100;

	double rho(std::vector<double> state) {
		return Mdot * exp(-2 * state[0] - state[1]) / (4 * PI * pow(r0, 2) * cs0);
	};
	double p(std::vector<double> state) {
		return rho(state) * T0 * exp(state[2]) / mb;
	};
	double cs(std::vector<double> state) {
		return sqrt(gamma * T0 * exp(state[2]) / mb);
	};
	double s(std::vector<double> state) {
		return log(T0 * exp(state[2] * pow(rho(state), 1 - gamma))) / (gamma - 1);
	};
	double dsdrho(std::vector<double> state) {
		return -1 / rho(state);
	};
	double dsdT(std::vector<double> state) {
		return 1 / ((1 - gamma) * T0 * exp(state[2]));
	};
	double dpds(std::vector<double> state) {
		return (gamma - 1) * p(state);
	};

};

int main() {
	double a = 10.;
	double b = 200;
	double beta = 1.;
	double g = 1.01;
	std::string fileName = "gammaheating.txt";
	std::ofstream output(fileName);
	
	//output.close();

	auto zeroFunc = [a,b,beta,g](double v0) {
		//std::cout << "v0=" << v0 << std::endl;
		GammaHeatingEquations eqs(a, b, beta, g, v0, true);
		Integrator solver(&eqs);
		double z = eqs.ZeroFunction(solver.getTLast(), solver.getDataLast());
		return z;
	};
	std::cout << "finding v0crit" << std::endl;
	clock_t timereq = clock();
	double v0crit = Find1DRoot(1.e-8, .9, zeroFunc, 1.e-8);
	timereq = clock() - timereq;
	std::cout << "v0crit=" << v0crit << ", time required: " << (float)timereq / CLOCKS_PER_SEC << std::endl;
	std::cout << zeroFunc(v0crit) << std::endl;
	//double v0crit = .00069885;

	for (double v0 = v0crit/2.; v0 < 2*v0crit; v0 += v0crit/5.) {
		std::cout << "making data for v0=" << v0 << std::endl;
		GammaHeatingEquations eqs(a, b, beta, g, v0, true);
		Integrator solver(&eqs);
		solver.writeToFile(&output);
	};

	std::cout << "making data for v0crit" << std::endl;
	GammaHeatingEquations eqs(a, b, beta, g, v0crit, true);
	Integrator solver(&eqs);
	solver.writeToFile(&output);
	output.close();
}
