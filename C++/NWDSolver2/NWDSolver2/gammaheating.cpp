#include <iostream> 
#include <fstream> 
#include <stdexcept>
#include <vector> 
#include <math.h>
#include <algorithm>
#include <ctime>

#include "integrator.h"
#include "rootfinder.h"
#include "EquationOfState.h"

class GammaHeatingEquations : public EoS {
public:
	GammaHeatingEquations(double a, double b, double beta, double gamma, double v0, bool critStop = false) : a(a), b(b), beta(beta), gamma(gamma), v0(v0), rmax(100.), critStop(critStop) {}
	int size() { return 3; }
	std::vector<double> getInitialState() { return { 0.0,log(v0),0.0 }; }

	double f1(std::vector<double> state) { 
		double val = 1. - exp(2. * state[1] - state[2]); 
		if (isinf(val)) { return copysign(LARGE_VALUE, val); };
		return val;
	};
	double f2(std::vector<double> state) { 
		double val = a * exp(-state[0] - state[2]) - 2. + exp(-state[0]-state[1]-state[2])*b*(1.-beta*exp(2*state[0]+6.*state[2]))*(gamma-1)/gamma;
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
		//std::cout << "differentials: dx=" << f1(state) << ", du=" << f2(state) << ", dw=" << dw << std::endl;
		if (critStop) {
			return { f1(state), f2(state), dw };
		};
		return { fabs(f1(state)), fabs(f2(state)), dw };
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
	bool critStop;
	const double LARGE_VALUE = 1.e100;
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
